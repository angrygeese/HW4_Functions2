# HW 4. Functions 
> *This is the repo for the third homework of the BI Python 2023 course*

# `protein_analyzer_tool.py`

## Title

This module contains the `protein_analyzer_tool.py` function that performs 5 operation on protein sequences. Operation are maintained using 5 secondary functions.


## Usage

To run the `protein_analyzer_tool.py`, first import it as module


```
import protein_analyzer_tool
```


and run its main function, `run_protein_analyzer_tool` function. This function provides interface for all 5 operations from `OPERATIONS` dictionary. Takes various number of positional arguments and one keyword-only argument:

- First `n` arguments - protein sequences;
- Latter positional argument - desired operation from list: "content_check", "seq_length", "protein_formula", "protein_mass", "charge";
- `abbrevition` keyword-only argument. Should be type integer, 1 for 1-letter abbreviation and 3 for 3-letter.

Returns tuple containing two list:

- `result` - list with operation results for each valid sequence;
- `corrupt_seqs` - list with non-valid sequences and their indices;


## Examples

Get molecular mass in g/mol for insulin:

```
>>> run_protein_analyzer_tool("AMALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN", "protein_mass", abbreviation=1)
All 1 sequence(s) processed successfully

(14090, [])
```


Get aminoacids content for various peptide:


```
>>> run_protein_analyzer_tool("AsnAspAspAsn", "content_check", abbreviation=3)
All 1 sequence(s) processed successfully

({'A': 0.0,
  'R': 0.0,
  'N': 50.0,
  'D': 50.0,
  'C': 0.0,
  'Q': 0.0,
  'E': 0.0,
  'G': 0.0,
  'H': 0.0,
  'I': 0.0,
  'L': 0.0,
  'K': 0.0,
  'M': 0.0,
  'F': 0.0,
  'P': 0.0,
  'S': 0.0,
  'T': 0.0,
  'W': 0.0,
  'Y': 0.0,
  'V': 0.0},
 [])
```


Such case provides no result:


```
>>> run_protein_analyzer_tool("", "protein_mass", abbreviation=1)
Processing result: [-]

0 sequence(s) out of 1 given have been processed successfully.
1 has been recognized as corrupted, i.e. non-protein

([], (0, ''))
```
Explanation of the latter is provided in next section


## Troubleshooting

`run_protein_analyzer_tool` raises errors in two cases:

*   Operation is not one from list: "content_check", "seq_length", "protein_formula", "protein_mass", "charge". If you are sure that input is correct, perform spell check.
*   Argument for `abrreviation` parameter is not integer from 1 or 3.


In other cases `run_protein_analyzer_tool` will not halt the execution. In other scenarios troubleshooting can be performed using second element in tuple returned by `run_protein_analyzer_tool`, `corrupt_seqs` list. This list contains sequences recognized as non-valid together with their indices in original sequence. in form of tuple `(<sequence_index>, <sequence>)`. Sequence is suggested to be non-valid in these cases:


*   If sequence is not type `str`. Other iterable objects are not supported by the time.
*   Sequence is empty string.

## Contacts
We hope our module provides useful tool for your work. If you encounter any errors, please mail one from our team: 

Belikova Angelina - kiit@gmail.com
Implemented: `protein_formula`, `protein_mass`, `seq_length`.

Aryuna Ayusheeva - aryuna.ayusheeva.1998@mail.ru
Implemented: ``aa_content_check`, `aa_chain_charge`.

Bredov Denis - d2707bredov@gmail.com
Teamlead. Implemented: `Mann_Whitney_U`, `decomposition`, `seq_transform`, `check_and_procees_seq`, `print_result`, `run_protein_analyzer_tool`.
