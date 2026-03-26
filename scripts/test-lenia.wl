(* --- Test Lenia script --- *)
PacletDirectoryLoad[FileNameJoin[{Directory[], "Lenia"}]];
Needs["Lenia`"];

Print["Running Lenia tests..."];
testReport = TestReport[FileNameJoin[{Directory[], "tests", "Lenia.wlt"}]];

If[testReport["AllTestsSucceeded"],
  Print["All ", testReport["TestsSucceededCount"], " tests passed!"],
  Print["Tests failed: ", testReport["TestsFailedCount"]];
  Scan[
    If[#["Outcome"] =!= "Success",
      Print["  FAIL: ", #["TestID"], " - ", #["Outcome"]]
    ] &,
    Values[testReport["TestResults"]]
  ];
  Exit[1]
];

(* --- Quick demo --- *)
Print["\nRunning a 64x64 demonstration..."];
gridSize = {64, 64};
R = 10;
(* Ring-shaped seed *)
grid = N @ Table[
  Module[{r = Sqrt[(i - 32)^2 + (j - 32)^2] / R},
    If[r < 1.0, Exp[-((r - 0.5)^2) / (2 * 0.15^2)], 0.0]
  ],
  {i, 1, 64}, {j, 1, 64}
];
result = Lenia[grid, 100, "Mu" -> 0.15, "Sigma" -> 0.015, "Radius" -> R];

If[Max[result] > 0.01,
  Print["Demo succeeded! Final max=", Max[result], " mean=", Mean[Flatten[result]]],
  Print["Demo failed: simulation decayed to zero."];
  Exit[1]
];

Print["Lenia implementation verified."];
