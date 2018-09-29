(* ::Package:: *)

(* ::Section:: *)
(*SMTrack*)


(* ::Subsection:: *)
(*Miscellaneous*)


BeginPackage["SMTrack`"];

Options[SMTrack] = {"segmented" -> False, "centroidW" -> 1.0, "sizeW" -> 0, "overlapW" -> 0, "subpixelLocalize" -> False,
"LoGkernel" -> 2, "threshold" -> 2.51, "constrainedSearch" -> False, "Dist" -> 15.0};


SMTrack::usage = "The package implements a robust single molecule tracking scheme. The procedure is somewhat similar to the
underlying algorithm proposed in Lineage Mapper (Chalfoun et al, Sci. Reports 2016). A novel implementation of a per-frame particle
jump distance has been incorporated, which relies on a mean distance(computed via a Delaunay Mesh) and its minimization such that a
maximum of 0 or 1 particle association is obtained between consecutive frames. The implementation is expected to work successfully
with numerous particle detection strategies. Another novel aspect is the inclusion of \"subpixel particle localization\" which can be
achieved using a two dimensional Gaussian-Fit";


(* mean separation between detected particles *)
meanParticleDist[centroids_]:= Module[{vertices,dist},
 vertices = MeshPrimitives[DelaunayMesh@centroids, 1]/.Line[x_]:> x;
 dist = Map[EuclideanDistance@@#&,vertices];
 {Mean@#,StandardDeviation@#}&@dist
];


(* where distDelaunay is the meanParticleDist@labeledMat *)
maxJumpDistance[centPrev_,centNew_,distDelaunay_]:= Module[{nearestFunc,rec,pts,\[ScriptCapitalA]},
 nearestFunc = Nearest@centPrev;
 rec = Flatten[Table[{i, Length@nearestFunc[#,{All, i}]},{i, 1, distDelaunay, 1}]&/@centNew, 1];
 pts = Cases[rec, {_,_?(# <= 1 &)}];
 (* \[ScriptCapitalA] = WeightedData[Keys@#,Exp[Values[#]/Total[Values@#]]]&@Counts[Part[pts,All,1]]; *)
 \[ScriptCapitalA] = WeightedData[Part[pts,All,1]];
 N@*Mean@\[ScriptCapitalA]
];


(*
segmentImage[image_Image,LoGKernelsize_,threshold_]:=MorphologicalComponents[
FillingTransform@MorphologicalBinarize[
ColorNegate@*ImageAdjust@LaplacianGaussianFilter[image,LoGKernelsize],
threshold]
];
*)


(*
segmentImage[image_Image,LoGKernelsize_]:=Module[{seg,seeds,bglabel},
seg =MorphologicalBinarize[
DeleteSmallComponents@Threshold[ColorNegate@LaplacianGaussianFilter[Blur[image],LoGKernelsize],
{"Hard","Cluster"}]
];
seeds=MaxDetect[DistanceTransform@seg];
seg=WatershedComponents[image,seeds,Method\[Rule]"Rainfall"];
bglabel=First@@MaximalBy[ComponentMeasurements[seg,"Area"],Values];
ArrayComponents[seg,2,bglabel\[Rule] 0]
];
*)


segmentImage[image_Image,LoGKernelsize_,threshold_]:=Module[{seg,seeds,bglabel},
 seg = MorphologicalBinarize[
  DeleteSmallComponents@Threshold[ColorNegate@Blur[LaplacianGaussianFilter[image,LoGKernelsize],threshold],
  {"Hard","Cluster"}]
 ];
 seeds = MaxDetect[DistanceTransform@seg];
 seg = WatershedComponents[image,seeds,Method->"Rainfall"];
 bglabel = First@@MaximalBy[ComponentMeasurements[seg,"Area"],Values];
 ArrayComponents[seg,2,bglabel-> 0]
];


modelFit[image_, mask_, shape_, box_] := Block[{pixelpos,pixelval,img,data,a,b,weights,data3D,
 mx,my,sx,sy,x,y,fm,dx,dy,cent,background,bestfit,intensityGuess,brightest,
 dimMask = Dimensions@mask},
 pixelpos = mask["NonzeroPositions"];
 pixelval = PixelValue[image, pixelpos];
 brightest = Max@pixelval;
 img = ReplacePixelValue[shape, Thread[PixelValuePositions[shape, 1] -> pixelval]];
 data = ImageData@ImagePad[img, 1];
 data3D = Flatten[MapIndexed[{First@#2, Last@#2, #1} &, data, {2}], 1];
 {dx,dy}= Dimensions@data;
 cent = N[(dx + dy)/2];
 intensityGuess = Max@data;
 fm = NonlinearModelFit[data3D,
  background+(a * Exp[-((x-mx)^2/(2*sx^2))-((y-my)^2/(2*sy^2))]),
 {{background,Automatic},{a,intensityGuess},{mx,cent/2},{my,cent/2},{sx,cent/4},{sy,cent/4}},
 {x, y}];
 bestfit = fm["BestFitParameters"];
 (* if any param is less than 0 then we run a constrained fit *)
 If[Length[Position[bestfit[[All,2]],x_/;x<0]]>0,
  fm = NonlinearModelFit[data3D,
 {background + (a * Exp[-((x-mx)^2/(2*sx^2)) - ((y-my)^2/(2*sy^2))]),
 {background > 0, a > 0, mx > 0, my > 0, sx > 0, sy > 0}},
 {{background,Automatic},{a,intensityGuess},{mx,cent/2},{my,cent/2},{sx,cent/4},{sy,cent/4}},
 {x,y}
 ];
];
 bestfit = fm["BestFitParameters"];
 If[brightest>0.9,
 weights =  data3D[[All,3]];
 weights = weights/.{x_/; x<0.90 -> 1.0,x_/;x<1.0 -> 0.0};
 (* set all intensities > 0.90 to zero weights *)
 fm = NonlinearModelFit[data3D,
 {background + (a*Exp[-((x-mx)^2/(2*sx^2)) - ((y-my)^2/(2*sy^2))]),
 {background > 0, a > 0 ,mx > 0, my > 0, sx > 0, sy > 0}},
 {{background,Automatic},{a,intensityGuess},{mx,cent/2},{my,cent/2},{sx,cent/4},{sy,cent/4}},
 {x,y}, Weights->weights
 ];
];
 bestfit = fm["BestFitParameters"];
 {background,a,mx,my,sx,sy} = {background,a,mx,my,sx,sy} /. fm["BestFitParameters"];
 (* {(Mean/@Transpose@box + ({mx,my} -{dx,dy}/2.0)),background,a,mx,my,sx,sy} *)
 If[AnyTrue[#,(#<0||#>First@dimMask||#>Last@dimMask)&],
  {Values@ComponentMeasurements[mask,"Centroid"],0,0,0,0,0,0}~FlattenAt~1,
  {#,background,a,mx,my,sx,sy}]&[(Mean/@Transpose@box+({mx,my}-{dx,dy}/2.0))]
];


funcGenerator[OptionsPattern[SMTrack]]:= Switch[OptionValue@"subpixelLocalize", True,
 (detectParticle[image_Image,LoGkernel_,thresh_]:= Block[{segImage, masks,boundingboxes,shapes},
 segImage = segmentImage[image,LoGkernel,thresh];
 {masks, shapes, boundingboxes} = Values@ComponentMeasurements[segImage,{"Mask","Shape","BoundingBox"}]\[Transpose];
 MapThread[modelFit[image, ##]&,{masks,shapes,boundingboxes}]
 ];
 (* prev and curr are props obtained from subpixelLocalization *)
 stackCorrespondence[prev_, curr_, opt:OptionsPattern[SMTrack]]:= Module[{costmat},
 costmat = costMatrix[prev[[All,1]],curr[[All,1]],opt];
 assignmentSubPixel[curr,costmat]
 ]), _ ,
 (detectParticle[image_Image,LoGkernel_,thresh_]:= segmentImage[image,LoGkernel,thresh];
 (* prev and curr are labeled matrices *)
 stackCorrespondence[prev_, curr_, opt: OptionsPattern[SMTrack]]:= Module[{costmat, currentMat,truelabels},
 costmat = costMatrix[prev,curr,opt];
 truelabels = Keys@ComponentMeasurements[prev,"Label"];
 assignmentLabelMat[curr,costmat,truelabels]
 ])
];


(* ::Subsection:: *)
(*Matrix Minimization*)


(* row wise minimums to determine which target cells are mapped from the source cells *)
rowwiseMins[costMat_]:= With[{constInfArray = ConstantArray[\[Infinity],Last@Dimensions[costMat]]},
 Rule@@@SparseArray[Unitize@Map[If[Min[#] == \[Infinity], constInfArray, # - Min@#]&,
  costMat], Automatic, 1]["NonzeroPositions"]
];


(* column wise minimums to determine cell mappings from current frame to the previous frame *)
colwiseMins[costMat_,groupingmetric_:(Last->First)]:= With[{constInfArray = ConstantArray[\[Infinity],First@Dimensions[costMat]]},
 GroupBy[
  SparseArray[Unitize@Map[If[Min[#] == \[Infinity], constInfArray ,#-Min[#]]&,
  costMat\[Transpose]], Automatic, 1]["NonzeroPositions"],
 groupingmetric, #]&
];


(* ::Subsection:: *)
(*Cost Matrix*)


overlapMatrix[seg1_,seg2_]:= Block[{keys1,keys2,mask,map,rules1,rules2},
 keys1 = Keys@ComponentMeasurements[seg1, "Label"];
 keys2 = Keys@ComponentMeasurements[seg2, "Label"];
 mask = Unitize[seg1*seg2];
 map = Normal@Counts@Thread[{SparseArray[mask*seg1]["NonzeroValues"],
 SparseArray[mask*seg2]["NonzeroValues"]}];
 {rules1, rules2}= Dispatch@Thread[# -> Range@Length@#]&/@{keys1, keys2};
 map[[All,1,1]] = map[[All,1,1]]/.rules1;
 map[[All,1,2]] = map[[All,1,2]]/.rules2;
 Normal@SparseArray[map,{Length@keys1,Length@keys2}]
];


(* for determining the overlapTerms *)
overlapCompiled = Compile[{{overlapmat, _Integer, 2},{prevList, _Real, 1}, {currList, _Real, 1}},
 1 - (overlapmat/(2.0 prevList) + (overlapmat\[Transpose] /(2.0 currList))\[Transpose]),
 CompilationTarget -> "C"
];


(* for determining the sizeTerm *)
sizeCompiled[prevList_, currList_] := With[{curr = currList},
 Abs[# - curr]/Clip[curr,{#,\[Infinity]}]&/@prevList
];


(* for determining the centroidTerm *)
centCompiled = Compile[{{centDiffMat, _Real, 2},{threshold, _Real}},
 Map[If[# >= threshold, 1., #/threshold]&, centDiffMat,{2}],
 CompilationTarget-> "C"
];


(* as the name implies, costMatrix generates the cost of traversing between an object @t and @t+1. More metrics e.g.
texture metrics can be added. In fact an arbitrary # of user defined metrics can be incorporated to compute the cost *)

costMatrix[Prev_,Curr_, OptionsPattern[SMTrack]]:= Module[{centroidPrev,centroidCurr,areaPrev,centroidDiffMat,areaCurr,nCol,
nRow,pos,mask,overlapMat,centroidW = OptionValue@"centroidW",spArraycentDiff,centroidTerm,maxCentDist,overlapTerm,sizeTerm,
spArrayOverlap,sizeW=OptionValue@"sizeW",overlapW=OptionValue@"overlapW",subpix=OptionValue@"subpixelLocalize"},

 If[subpix,
  {centroidPrev,centroidCurr} = {Prev,Curr},
  {centroidPrev,areaPrev}=Values@ComponentMeasurements[Prev,{"Centroid","Area"}]\[Transpose];
  {centroidCurr,areaCurr}=Values@ComponentMeasurements[Curr,{"Centroid","Area"}]\[Transpose];
  ];
  {nRow,nCol} = Length/@{centroidPrev,centroidCurr};
  maxCentDist = If[OptionValue["constrainedSearch"], OptionValue["Dist"],
    maxJumpDistance[#,centroidCurr,First@meanParticleDist[#]]&@centroidPrev
   ];
  centroidDiffMat = DistanceMatrix[N@centroidPrev,N@centroidCurr];
  centroidTerm = centCompiled[centroidDiffMat,maxCentDist];
  spArraycentDiff = SparseArray[UnitStep[maxCentDist - centroidDiffMat], Automatic, 0];

  pos = If[overlapW > 0,(overlapMat = overlapMatrix[Prev,Curr];
  overlapTerm = overlapCompiled[overlapMat,areaPrev,areaCurr]; (* compute overlapTerm *)
  spArrayOverlap = SparseArray[overlapTerm,Automatic,1.];
  (* positions where overlaps occur or centroids within maxCentDist *)
  spArrayOverlap["NonzeroPositions"]~Union~spArraycentDiff["NonzeroPositions"]),
   (overlapTerm = 0; spArraycentDiff["NonzeroPositions"])
  ];
  (* compute sizeTerm *)
  sizeTerm = If[sizeW > 0, areaPrev~sizeCompiled~areaCurr, 0];
  (* creating the mask for costMat *)
  mask = SparseArray[pos -> 1, {nRow,nCol}, \[Infinity]];
  mask*(overlapW*overlapTerm + centroidW*centroidTerm + sizeW*sizeTerm)
];


(* ::Subsection:: *)
(*Assignment Problem*)


assignmentHelper[costMat_,truePrevKeys_]:=Module[{rmins,otherassoc,rules,assignmentList,indices,
artificialInds,previnds,realindices,ruleAssigned,currentassigned},
  rmins = rowwiseMins@costMat;
  (* other possible associations using columnwise mins *)
  otherassoc = Sort@Flatten@KeyValueMap[Thread@*Rule,
  DeleteCases[colwiseMins[costMat]@Identity,x_/;Length@x>1]];
  rules = Sort@DeleteDuplicates[rmins~Join~otherassoc];
  artificialInds = rules/.Rule->List;
  previnds = Part[ truePrevKeys,artificialInds[[All,1]] ];
  realindices = Transpose[{previnds,Part[artificialInds,All,2]}];(*{true label prev, label current}*)
  (* create the graph with initialized weights for the assignment problem *)
  assignmentList = Block[{p,c,edges,edgeweights,graph,assignments},
  edges = Subscript[p,First@#]->Subscript[c,Last@#]&/@realindices;
  edgeweights = costMat[[Sequence@@#]]&/@artificialInds;
  graph = Graph[edges,EdgeWeight->edgeweights];
  assignments = FindIndependentEdgeSet@graph;
  Replace[assignments, HoldPattern[Subscript[p,x_]\[DirectedEdge] Subscript[c,y_]]:>{x,y},{1}]
  ];
  {artificialInds[[All,1]],assignmentList}
 ];


assignmentLabelMat[segCurr_,costMat_,truePrevKeys_]:= Module[{segmentCurr= segCurr, currframelabels,
 newlabels, ruleAssigned, currentassigned,currentunassigned, maxlabelprev, newcellAssignmentRules,
 allAssignmentRules,assignmentsList,artificialInds},

 currframelabels = Keys@ComponentMeasurements[segCurr,"Label"];
 {artificialInds,assignmentsList} = assignmentHelper[costMat,truePrevKeys];
 ruleAssigned = Reverse[Rule@@@assignmentsList,{2}];
 currentassigned = Part[assignmentsList, All, 2];
 currentunassigned = Complement[currframelabels,currentassigned]; (* new spots *)
 maxlabelprev = Max@truePrevKeys;
 newlabels = Range[maxlabelprev+1,maxlabelprev+Length@currentunassigned];
 newcellAssignmentRules = Thread[currentunassigned-> newlabels];
 allAssignmentRules = Dispatch@SortBy[newcellAssignmentRules~Join~ruleAssigned, First];
 Replace[segCurr,allAssignmentRules,{2}]
];


DownValues[caten]={
 HoldPattern[caten[{},{}]]:>{},
 HoldPattern[caten[list1:{{__?NumberQ},__},list2:{{__?NumberQ},__}]]:>Join[{list1},{list2}],
 HoldPattern[caten[list1_,list2:{{__?NumberQ},__}]]:>Join[list1,{list2}],
 HoldPattern[caten[list1:{{__?NumberQ},__},list2_]]:>Join[{list1},list2],
 HoldPattern[caten[list1:{{{__?NumberQ},__}..},list2_]]:>Join[list1,list2],
 HoldPattern[caten[list1_,list2:{{{__?NumberQ},__}..}]]:>Join[list1,list2]
};


assignmentSubPixel[Curr_,costMat_]:= Module[{ newlabels,assignmentsList,currentassigned,currentunassigned,
 maxlabelprev,rules, dim = Dimensions@costMat,parent,artificialInds},

 (* add tracked label points to seeds *)
 {artificialInds,assignmentsList} = assignmentHelper[costMat,currentKeyVector];
 assignmentsList = SortBy[assignmentsList,First];
 {parent,currentassigned } = {#[[All,1]],#[[All,2]]}&[assignmentsList];

 MapThread[(seeds[#1] = Join[seeds[#1],{#2}])&,{parent,Part[Curr,currentassigned]}];

 (* new spots added to seeds and currentKeyVector *)
 currentunassigned = Complement[Range[Last@dim],currentassigned];
 maxlabelprev = Max@currentKeyVector;
 newlabels = Range[maxlabelprev + 1, maxlabelprev + Length@currentunassigned];

 If[newlabels != {},
  AssociateTo[seeds,
   MapAt[List,
    Thread[newlabels -> Part[Curr,currentunassigned]],{All,2}]
  ]
];

 (* remove untracked parents from currentKeyVector *)
 currentKeyVector = DeleteCases[currentKeyVector,Alternatives@@(currentKeyVector~Complement~parent)];
 currentKeyVector = Join[currentKeyVector,newlabels];

 If[currentassigned != {},
  (Curr[[currentassigned]]~caten~Curr[[currentunassigned]]),
   Curr[[currentunassigned]]
  ]
];


(* ::Subsection:: *)
(*Main Function*)
SMTrack[filename_, opt: OptionsPattern[]]:= Module[{segmented = OptionValue["segmented"], input,
 subpixloc = OptionValue@"subpixelLocalize", imports= Import@filename, logKernel = OptionValue@"LoGkernel",
 thresh = OptionValue@"threshold"},

 funcGenerator["subpixelLocalize" -> subpixloc];

 input = Switch[segmented, False,
 If[subpixloc == False,
   ParallelTable[Quiet@detectParticle[i,logKernel,thresh],{i,imports}]
  ], _ ,imports];

If[subpixloc,
 currentKeyVector = Range[Length@*First@input];
 seeds = <|Thread[currentKeyVector -> Partition[First@input,1]]|>
 ];

 FoldList[stackCorrespondence[##,opt]&,First@input,Rest@input]
];


(* ::Subsection:: *)
(*End Package [ ]*)


Begin["`Private`"];
currentKeyVector;
seeds;
End[];


EndPackage[];
