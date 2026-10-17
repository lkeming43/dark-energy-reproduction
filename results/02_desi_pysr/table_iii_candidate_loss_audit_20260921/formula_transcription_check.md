# Table III formula-transcription check

## Source check

The official arXiv v2 source archive was inspected directly. The unusual notation is present in the authors' LaTeX source; it is not introduced by OCR, the screenshot, or the arXiv HTML conversion.

Relevant source rows include:

- complexity 8: `...^{2} - -3.47`
- complexity 13: `... + 6.56 - - \frac{0.02}{z-0.76}`
- complexity 15: `\left(3.68 - z\times\right)^{2.37}`
- complexity 17: `... + 6.56 - - \frac{0.0169}{...}`

## Interpretation used in the audit

The two consecutive minus signs were interpreted algebraically:

`a - (-b) = a + b`.

This is consistent with raw symbolic-regression output containing negative subexpressions. It is also strongly favored by a direct sign test:

| row | double-minus interpreted as plus | interpreted as a single minus | paper loss |
| --- | ---: | ---: | ---: |
| complexity 8 | 0.248055 | 5455.896400 | 0.24 |
| complexity 13 | 1.299872 | 2.696304 | 0.02 |
| complexity 17 | 0.054485 | 0.961693 | 0.013 |

For complexity 8, the plus interpretation nearly reproduces the reported loss, whereas treating the duplicated sign as a typographical single minus fails catastrophically. The same interpretation also gives the closer result for complexities 13 and 17.

## Complexity 15

Complexity 15 is genuinely malformed because the multiplication sign after `z` has no right-hand operand. The audit removes only that dangling multiplication sign and labels the result as a minimal repair:

`(3.68-z)^2.37 + 6.56 + 0.0169/[0.862*(z-0.763)]`.

This repair is supported by complexity 16, whose first term

`21.9*(1-0.272*z)^2.37`

is a rounded algebraic rewriting of `(3.68-z)^2.37`, since `3.68^2.37` is approximately `21.9` and `1/3.68` is approximately `0.272`. Nevertheless, the missing source token means complexity 15 cannot be treated as an exact reproduction.

## Impact on the main conclusion

Formula parsing does not explain the principal loss discrepancy:

- complexities 1--5 are syntactically clear;
- the selected complexity-5 expression `28.42/(1.68^z)` is unambiguous and equivalent to an exponential;
- its recomputed public-covariance loss is 0.664026 rather than 0.30;
- even refitting the two exponential coefficients gives a minimum loss of about 0.407375 with the public covariance, still above 0.30;
- complexity 9 is parsed literally and does reproduce the displayed paper loss.

Therefore one malformed high-complexity row should be excluded, but it does not invalidate the earlier conclusion that the Table III loss column is not uniformly reproducible from the published formulas, public data, and stated Eq. (9).

