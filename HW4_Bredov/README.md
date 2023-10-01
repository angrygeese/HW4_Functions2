# HW 4. Functions 
> *This is the repo for the third homework of the BI Python 2023 course*

# `protein_analyzer_tool.py`

## Title

This module contains the `protein_analyzer_tool.py` function that performs 5 basic operations on protein sequences. `protein_analyzer_tool` uses 5 secondary functions.

## Overview
- `protein_mass` - returns molecular weight of protein in g/mol. The function takes string and returns integer.
- `protein_formula` - returns the molecular formula of protein. The function takes string and returns dictionary with data of the protein atomic composition. The output contains the number of C, H, N, O, and S (optionally) atoms in sequence.
## Usage

To run the `protein_analyzer_tool.py`, first import it as module


```
import protein_analyzer_tool
```


and its main function, `run_protein_analyzer_tool` function. This function takes 3 arguments:



## Options
## Examples

```
>>> run_protein_analyzer_tool("AMALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN", "protein_mass", abbreviation=1)
All 1 sequence(s) processed successfully

(14090, [])
```

```
>>> run_protein_analyzer_tool("", "protein_mass", abbreviation=1)
Processing result: [-]

0 sequence(s) out of 1 given have been processed successfully.
1 has been recognized as corrupted, i.e. non-protein

([], (0, ''))
```
Further description of this case is given in next section


## Troubleshooting

`run_protein_analyzer_tool` raises errors in two cases:

*   Operation is not one from list: "content_check", "seq_length", "protein_formula", "protein_mass", "charge". If you are sure that input is correct, perform spell check.
*   Argument for `abrreviation` parameter is not integer from 1 or 3.


In other cases `run_protein_analyzer_tool` will not halt the execution. In other scenarios troubleshooting can be performed using second element in tuple returned by `run_protein_analyzer_tool`, `corrupt_seqs` list. This list contains sequences recognized as non-valid together with their indices in original sequence. in form of tuple `(<sequence_index>, <sequence>)`. Sequence is suggested to be non-valid in these cases:


*   If sequence is not type `str`. Other iterable objects are not supported by the time.
*   Sequence is empty string.

## Contacts
Bredov Denis - d2707bredov@gmail.com
Belikova Angelina - kiit@gmail.com
