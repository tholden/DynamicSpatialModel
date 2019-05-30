DynareTransformationEngine
==========================

A set of mod files for automatically transforming variables into logs or logits as appropriate, and for adding similarly transformed AR(1) shock processes.

Compare `example_original.mod` (not using the engine) and `example.mod` (using the engine) to see basic usage.

For leads greater than 1, insert: `@#define MaximumLead = 2` at the top of your MOD file and then use `..._LEAD2`, `..._LEAD3`, etc.

For lags greater than 1, insert: `@#define MaximumLag = 2` at the top of your MOD file and then use `...LAG2`, `...LAG3`, etc.
