function controlDefs = define_control()

%   "Variables",            "VarMin",           "VarMax",   "Conditions"                "Transformations"
controlDefs = {...
    "ControlSpan",         [NaN,NaN],          [NaN,NaN]    "> Previous",               "~";
    "ControlChord",         NaN,                NaN,        "~",                        "~"};