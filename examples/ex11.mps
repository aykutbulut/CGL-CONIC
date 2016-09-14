NAME
* Example from Goez's dissertation, section 5.2.
OBJSENSE
    MIN
ROWS
 N  OBJ
 E  C1
 E  C2
 E  C3
 E  C4
 E  C5
COLUMNS
    X1        OBJ       3e+0
    X1        C1        9e+0
    X1        C2        1e+0
    X1        C3        1e+0
    X1        C4        1e+0
    0         'MARKER'                 'INTORG'
    X2        OBJ       1e+0
    X2        C1        1e+0
    X2        C4        1e+0
    X2        C5        1e+0
    1         'MARKER'                 'INTEND'
    X3        OBJ       2e+0
    X3        C1        1e+0
    X3        C3        1e+0
    X4        OBJ       1e+0
    X4        C1        1e+0
    X5        OBJ       0
    X5        C1        1e+0
    X5        C2        9e+0
    X5        C3        1e+0
    X5        C4        1e+0
    X5        C5        1e+0
    X6        OBJ       0
    X6        C2        1e+0
    X6        C4        1e+0
    X6        C5        1e+0
    2         'MARKER'                 'INTORG'
    X7        OBJ       1e+0
    X7        C2        1e+0
    X7        C4        1e+0
    3         'MARKER'                 'INTEND'
    X8        OBJ       1e+0
    X8        C2        1e+0
    X8        C3        1e+0
    X8        C5        1e+0
RHS
    rhs       C1        1e+1
    rhs       C2        1e+1
    rhs       C3        2e+0
    rhs       C4        1e+0
    rhs       C5        1e+0
RANGES
BOUNDS
 FR bound     X1
 FR bound     X2
 FR bound     X3
 FR bound     X4
 FR bound     X5
 FR bound     X6
 FR bound     X7
 FR bound     X8
CSECTION      K1        0              QUAD
    X1
    X2
    X3
    X4
CSECTION      K2        0              QUAD
    X5
    X6
    X7
    X8
ENDATA
