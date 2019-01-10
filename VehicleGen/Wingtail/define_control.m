function controlDefs = define_control()

controlDefs = {...
    "Variables",        "VarMin",       "VarMax",       "Conditions",	"Transformations",  "Optimise/Hold";
    "ControlSpan",     [NaN,NaN],       [NaN,NaN],      "> Previous",	"~",                "Hold";
    "ControlChord",     NaN,             NaN,           "~",            "~",                "Hold"};

controlDefs = optorhold(controlDefs);