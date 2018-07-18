(* ::Package:: *)

(* ::Section:: *)
(*Spot Tracking*)


BeginPackage["SMTrack`", {"segmentParticles`"}];


SMTrack::usage = "implements a robust single molecule tracking algorithm. The procedure is similar to the underlying algorithm proposed in Lineage Mapper (Chalfoun et al). However,
a Kalman predictor has been incorporated into the tracking scheme to address the case of occluding particles when the motion is linear. The implementation is expected to work successfully
with numerous particle detection strategies.";


(* ::Subsection:: *)
(*misc*)


Options[SMTrack] = {"segmented" -> False, "centroidW" -> 1.0, "maxCentDist" -> 10.0};


(* mean separation between detected particles *)
particleSeparation[labeledMat_]:= Module[{vertices,dist},
vertices = MeshPrimitives[DelaunayMesh@Values@ComponentMeasurements[labeledMat,"Centroid"],1]/.Line[x_]:> x;
dist = Map[EuclideanDistance@@#&,vertices];
{Mean@#,StandardDeviation@#}&@dist
];


(* ::Subsection:: *)
(*cost matrix*)


(* for determining the centroidTerm *)
centCompiled = Compile[{{centDiffMat, _Real, 2},{threshold, _Real}}, 
Map[If[# >= threshold, 1., #/threshold]&, centDiffMat,{2}],
CompilationTarget-> "C"
];


(* as the name implies, costMatrix generates the cost of traversing between an object @t and @t+1. More metrics e.g. texture metrics can be added.
In fact an arbitrary # of user defined metrics can be incorporated to compute the cost *)
costMatrix[segPrev_,segCurr_, OptionsPattern[SMTrack]]:= Module[{centroidPrev,centroidCurr, centroidDiffMat,
areaPrev,areaCurr, nCol,nRow,pos,centroidTerm,spArraycentDiff,mask, centroidW = OptionValue@"centroidW", 
maxCentDist = OptionValue@"maxCentDist"},

{centroidPrev,areaPrev} = Values@ComponentMeasurements[segPrev,{"Centroid","Area"}]\[Transpose];
{centroidCurr,areaCurr} = Values@ComponentMeasurements[segCurr,{"Centroid","Area"}]\[Transpose];
{nRow,nCol} = {Length@areaPrev,Length@areaCurr};

centroidDiffMat = DistanceMatrix[N@centroidPrev,N@centroidCurr]; (* calculating pairwise distance between centroids in frames t, t+1 *)
centroidTerm = centCompiled[centroidDiffMat,maxCentDist]; (* optimized to compute centroidTerm *)

spArraycentDiff = SparseArray[UnitStep[maxCentDist - centroidDiffMat], Automatic, 0]; (* finding positions in centDiffMat for dist within maxCentDist *)
pos = spArraycentDiff["NonzeroPositions"]; (* positions where overlaps occur or centroids within maxCentDist *)
mask = SparseArray[pos -> 1,{nRow,nCol},\[Infinity]]; (* creating the mask for costMat *)
mask*(centroidW*centroidTerm)
];


(* ::Subsection:: *)
(*matrix minimization*)


(* row wise minimums to determine which target cells are mapped from the source cells *)
rowwiseMins[costMat_]:= With[{constInfArray = ConstantArray[\[Infinity],Last@Dimensions[costMat]]},
Rule@@@SparseArray[Unitize@Map[If[Min[#] == \[Infinity], constInfArray, # - Min@#]&,costMat], Automatic, 1]["NonzeroPositions"]
];


(* column wise minimums to determine cell mappings from current frame to the previous frame *)
colwiseMins[costMat_,groupingmetric_:(Last->First)]:= With[{constInfArray = ConstantArray[\[Infinity],First@Dimensions[costMat]]},
GroupBy[
SparseArray[Unitize@Map[If[Min[#] == \[Infinity],constInfArray ,# - Min[#]]& , costMat\[Transpose]], Automatic, 1]["NonzeroPositions"],
groupingmetric, #]&
];


(* ::Subsection:: *)
(*assignment problem*)


(* graph theoretic approach to finding assignments *)
assignmentFunc[segCurr_,costMat_,truePrevKeys_]:= Module[{segmentCurr= segCurr, currframelabels, realindices,
 artificialInds, newlabels, assignmentsList, ruleAssigned, currentassigned,currentunassigned, maxlabelprev, newcellAssignmentRules,
 allAssignmentRules, previnds, otherassoc, rules, rmins},
currframelabels = Keys@ComponentMeasurements[segCurr,"Label"];
rmins = rowwiseMins@costMat;
(* other possible associations using columnwise mins *)
otherassoc = Sort@Flatten@KeyValueMap[Thread@*Rule, DeleteCases[colwiseMins[costMat]@Identity,x_/;Length@x>1]];
rules = Sort@DeleteDuplicates[Join[rmins,otherassoc]]; 
artificialInds = rules/.Rule -> List;
previnds = Part[ truePrevKeys,artificialInds[[All,1]] ];
realindices = Transpose[{previnds,Part[artificialInds,All,2]}];

(* create the graph with initialized weights for the assignment problem *)
assignmentsList = Block[{p,c,edges,edgeweights,graph,assignments},
edges=Subscript[p, #[[1]]]-> Subscript[c, #[[2]]]&/@realindices;
edgeweights = costMat[[Sequence@@#]]&/@artificialInds;
graph = Graph[edges,EdgeWeight->edgeweights];
assignments= FindIndependentEdgeSet@graph;
Replace[assignments, HoldPattern[Subscript[p, x_] \[DirectedEdge] Subscript[c, y_]]:>{x,y},{1}]
];

ruleAssigned = Reverse[Rule@@@assignmentsList,{2}];
currentassigned = Part[assignmentsList,All,2];
currentunassigned = Complement[currframelabels,currentassigned]; (* these are cells that arrived into FOV *)
maxlabelprev = Max@truePrevKeys;
newlabels = Range[maxlabelprev+1,maxlabelprev+Length@currentunassigned];
newcellAssignmentRules = Thread[currentunassigned-> newlabels];
allAssignmentRules = Dispatch@SortBy[newcellAssignmentRules~Join~ruleAssigned,First];
Replace[segCurr,allAssignmentRules,{2}] 
]


(* ::Subsection:: *)
(*helper( ) and main( )*)


(* Helper Function to Main *)
stackCorrespondence[prev_, curr_]:=Module[{costmat, currentMat,truelabels},
(* ensure the current stack is segmented into a labeled matrix *)
currentMat = If[Head@curr === Image, segmentImage[curr], curr];
(* compute the cost matrix *)
costmat = costMatrix[prev,currentMat];
(* assignments and replacements *)
truelabels = Keys@ComponentMeasurements[prev,"Label"];
assignmentFunc[currentMat,costmat,truelabels]
];


(* Main *)
SMTrack[filename_, opt: OptionsPattern[]]:= Module[{option = OptionValue["segmented"],imports, Prev},
imports= Import@filename;
Prev = Switch[option, False, segmentImage@*First@imports, _, First@imports];
FoldList[stackCorrespondence[##,opt]&,Prev,Rest@imports ]
];


(* ::Subsection:: *)
(*End Package [ ]*)


Begin["`Private`"];
End[];


EndPackage[];
