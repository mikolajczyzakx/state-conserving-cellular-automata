# One dimensional state-conserving Cellular Automata with the four-cell neighborhood

This repository contains files with lookup tables (LUTs) for one-dimensional state-conserving cellular automata with **four-cell neighborhood** and **three states**.

---

## LUT Definition

A **LUT** (lookup table) is the complete set of all local transition rules of a cellular automaton — i.e. the value of the transition function `f` for every possible neighborhood configuration:

```
f(s1, s2, s3, s4) → {0, 1, 2},   where  s1, s2, s3, s4 ∈ {0, 1, 2}
```

For neighborhood size 4 and 3 states there are **3⁴ = 81** possible neighborhood configurations, so a LUT consists of **81 local rules**, listed in the order shown below:

```
f(2,2,2,2), f(2,2,2,1), ... , f(0,0,0,1), f(0,0,0,0)
```

Each line in a file is one complete LUT — a sequence of 81 output values in the order above.

### Numeric interpretation

This ordering is equivalent to writing the LUT as a number in **base 3** (ternary positional numeral system), where:

- the **most significant digit** is `f(2,2,2,2)`,
- the **least significant digit** is `f(0,0,0,0)`.

So a LUT can be read as an 81-digit number in base 3.

This gives every LUT a unique numeric identifier and a canonical total ordering.

---

## Files

| File | Description |
|---|---|
| [LUTs_all.txt](./LUTs/LUTs_all.txt) | All 6312 LUTs for the automata — complete set |
| [LUTs_type0.txt](./LUTs/LUTs_type0.txt) | 3155 Type 0 automata |
| [LUTs_swap.txt](./LUTs/LUTs_swap.txt) | 2348 Swap automata |
| [LUTs_bypass.txt](./LUTs/LUTs_bypass.txt) | 48 Bypass automata |
| [LUTs_mix.txt](./LUTs/LUTs_mix.txt) | 807 Swap+Bypass automata |