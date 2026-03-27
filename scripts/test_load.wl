PacletDirectoryLoad[FileNameJoin[{Directory[], "Lenia"}]];
Needs["Lenia`"];
Print["Calling initRustFunction..."];
Print[ "Library path: ", Lenia`Private`leniaRustLibPath[] ];
result = Lenia`Private`initRustFunction[];
Print["Result: ", result];
Print["Done."];
