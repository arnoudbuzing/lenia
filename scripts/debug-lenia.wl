(* Test LeniaSeed types *)
PacletDirectoryLoad[FileNameJoin[{Directory[], "Lenia"}]];
Needs["Lenia`"];

types = {"Ring", "MultiRing", "RandomOrganism", "Constellation", "Asymmetric"};

Do[
  grid = LeniaSeed[128, type];
  Print[type, ": grid max=", Max[grid], " mean=", Mean[Flatten[grid]]];
  result = Lenia[grid, 100, "Mu" -> 0.15, "Sigma" -> 0.015];
  Print["  After 100 steps: max=", Max[result], " mean=", Mean[Flatten[result]]],
  {type, types}
];

(* Test random selection *)
Print["\nRandom seed:"];
grid = LeniaSeed[128];
Print["  grid max=", Max[grid], " dims=", Dimensions[grid]];
result = Lenia[grid, 50];
Print["  After 50 steps: max=", Max[result], " mean=", Mean[Flatten[result]]];
