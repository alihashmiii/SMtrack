(* ::Package:: *)

(* ::Section:: *)
(*Spot Detection*)


BeginPackage["detectParticles`"]


(* ::Subsection:: *)
(*Segment Image*)


segmentImage[image_Image,LoGKernelsize_,threshold_]:=MorphologicalComponents[
FillingTransform@MorphologicalBinarize[ColorNegate@ImageAdjust@LaplacianGaussianFilter[image,LoGKernelsize],threshold]
];


Begin["`Private`"]
End[] 
EndPackage[]
