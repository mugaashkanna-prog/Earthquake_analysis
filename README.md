# Earthquake Magnitude Range Analysis Using a Wavelet Tree

A Python-only DSA lab case study. A hand-written **Wavelet Tree** answers
magnitude queries (k-th largest/smallest, threshold counts, median,
distribution, period comparison) on any date range in `O(log σ)` time,
without re-sorting or re-scanning the data.

## Files

| File | Purpose |
|---|---|
| `earthquake_analysis.py` | Complete program: loader, validation, compression, Wavelet Tree, queries, menu CLI, `--demo` mode |
| `test_earthquake_analysis.py` | 18 unit tests; every answer is compared with plain Python sort/count |
| `sample_earthquakes.csv` | Small hand-made test file (28 valid rows + 5 deliberately bad rows) |
| `generate_synthetic_data.py` | Makes a large **random** CSV (not real data) for speed testing |

Requires Python 3.8+ and only the standard library (`csv`, `datetime`, `bisect`, `time`).

## 1. Getting the real USGS dataset

The sample file is illustrative (hand-typed values loosely based on well-known
events, **not authoritative data**). For your submission, use real USGS data (≥10,000 rows):

**Option A – website**
1. Open <https://earthquake.usgs.gov/earthquakes/search/>
2. Set *Output Options* → **CSV**.
3. Set *Date & Time* (e.g. 2015-01-01 to 2024-12-31) and *Magnitude* minimum (e.g. 4.5; or 2.5 for a worldwide set).
4. Click **Search** and save the file, e.g. as `usgs_earthquakes.csv`.

**Option B – download link (same service)**
```bash
curl -o usgs_2015_2024.csv "https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv&starttime=2015-01-01&endtime=2024-12-31&minmagnitude=4.5"
```
USGS refuses a single request returning more than **20,000 events**. If you
see an error, download in pieces (one year, or one magnitude band, at a time)
and join them by keeping one header line:
```bash
head -n 1 part1.csv > usgs_all.csv
tail -n +2 -q part1.csv part2.csv part3.csv >> usgs_all.csv
```
The program reads these columns by name: `time, latitude, longitude, depth, mag, place, id`
(only `time` and `mag` are mandatory). The standard USGS CSV already has them.
No other preparation is needed – the program validates, de-duplicates and sorts it.

## 2. Running

```bash
python earthquake_analysis.py usgs_all.csv         # interactive menu
python earthquake_analysis.py                      # uses sample_earthquakes.csv
python earthquake_analysis.py --demo               # non-interactive demo + self-check
python test_earthquake_analysis.py                 # run all tests
python generate_synthetic_data.py 20000 big.csv    # optional large random test file
```
Dates are `YYYY-MM-DD`, or just `YYYY` (start = 1 Jan, end = 31 Dec). Both ends are **inclusive**.

## 3. Sample input / output

Sample file, menu options 1, 3, 4, 5, 6 (inputs shown after the `:`):

```
Loaded 'sample_earthquakes.csv' and built the Wavelet Tree in 0.005 s
Rows read: 33 | kept: 28 | bad magnitude: 3 | bad time: 1 | duplicate id: 1
Time span : 2019-03-14  ->  2022-12-20
Magnitudes: min 2.70, max 8.20, 22 distinct values (sigma)

Choose: 1                                  (Feature 1: k-th largest)
Start date: 2020   End date: 2022   k: 3
Period 2020-01-01 to 2022-12-31 -> 23 earthquake(s)
k = 3 -> k-th largest magnitude: 7.70

Choose: 3                                  (Feature 3: frequency)
Start 2019  End 2022  Threshold: 5.0
Earthquakes with magnitude >= 5.0: 20

Choose: 4                                  (Feature 4: median)
Start 2020  End 2020
Median magnitude: 5.400                    (10 quakes -> mean of 5.3 and 5.5)

Choose: 5                                  (Feature 5: distribution, year 2020)
  Below 3.0       : 1
  3.0 - 4.9       : 3
  5.0 - 6.9       : 3
  7.0 and above   : 3

Choose: 6                                  (Feature 6: A = 2019-2020, B = 2021-2022, threshold 7.0)
                             A: 2019-01-01..2020-12-31 B: 2021-01-01..2022-12-31
Number of earthquakes        15                     13
Maximum magnitude            7.70                   8.20
Median magnitude             5.20                   7.20
Count with mag >= 7.0        4                      7
```

Error handling (real output):
```
Invalid input: day is out of range for month. Please try again.     <- date 2030-02-31
Invalid k: it must be between 1 and 10.                             <- k = 99
No earthquakes in this period - nothing to rank.                     <- empty range
No earthquakes in this period - median is undefined.
```

`--demo` output on the 20,000-row synthetic file (each answer cross-checked):
```
Period 2020-01-01 .. 2020-12-31  (2004 quakes, index range [10098, 12102))
   3rd largest      wavelet=6.5  plain-python=6.5  [OK]
   count >= 5.0     wavelet=44   plain-python=44   [OK]
   median           wavelet=2.9  plain-python=2.9  [OK]
```

## 4. The Wavelet Tree explained

**Problem it solves.** Chronological sorting turns "a date range" into a
contiguous slice `[left, right)` of the magnitude array. We then need order
statistics (k-th smallest, rank, median) of that *slice*, for many different
slices. A Wavelet Tree answers these directly.

**Step 0 – coordinate compression.** Magnitudes are floats, so we replace each
by its rank among the sorted *distinct* values:
`[4.5, 7.1, 4.5, 3.0, 5.0]` → distinct `[3.0, 4.5, 5.0, 7.1]` → ranks `[1, 3, 1, 0, 2]`.
Equal magnitudes share a rank (duplicates work) and rank order = magnitude order.
Now values are small integers `0..σ-1`.

**Construction (worked example).** Sequence `[1,3,1,0,2]`, rank range `[0,3]`, `mid = 1`.

```
Root  range [0,3]   sequence  1 3 1 0 2
      goes left?              L R L L R        (rank <= 1 -> left)
      left_prefix = [0, 1, 1, 2, 3, 3]
        ├── Left  [0,1]  seq 1 1 0   mid=0   left_prefix=[0,0,0,1]
        │     ├── leaf 0   (rank 0 = 3.0)
        │     └── leaf 1   (rank 1 = 4.5)
        └── Right [2,3]  seq 3 2      mid=2   left_prefix=[0,0,1]
              ├── leaf 2   (rank 2 = 5.0)
              └── leaf 3   (rank 3 = 7.1)
```
Each node stores only `left_prefix[i]` = "how many of my first `i` items went left".
Children's sequences are never stored – they are implied by the arrays above.

**Mapping a range downward.** A range `[l, r)` at a node becomes
- left child: `[left_prefix[l], left_prefix[r])`
- right child: `[l - left_prefix[l], r - left_prefix[r])`

**k-th smallest, worked example:** 2nd smallest in positions `[1,5)` (ranks `3,1,0,2`).
Root: `left_prefix[1]=1, left_prefix[5]=3` → 2 items go left; `k=2 ≤ 2` → go left with range `[1,3)`.
Left child: `left_prefix[1]=0, left_prefix[3]=1` → 1 item goes left; `k=2 > 1` → `k=1`, go right with range `[1,2)`.
Leaf with rank 1 reached → magnitude **4.5**. ✔ (the slice sorted is 3.0, 4.5, 5.0, 7.1)

**Other queries**
- *k-th largest* = `(length − k + 1)`-th smallest.
- *count of magnitude ≥ T*: `bisect` finds the first rank with value ≥ T, then
  `count_less_than_rank` walks down once (adding the "went left" count whenever it
  must go right); answer = `length − count_less`.
- *Distribution*: boundaries 3.0, 5.0, 7.0 with **half-open** buckets
  `<3.0`, `[3.0,5.0)`, `[5.0,7.0)`, `≥7.0` – every record is in exactly one bucket
  (so 4.95 is in "3.0–4.9", 5.0 in "5.0–6.9"). The four counts always sum to the range size.
- *Median*: odd count → middle; even → average of the two middle k-th smallests.
- *Dates → indices*: `bisect_left(timestamps, start)` and
  `bisect_left(timestamps, end + 1 day)` give `[left, right)`.

## 5. Complexity

n = records, σ = distinct magnitudes, depth ≈ ⌈log₂ σ⌉.

| Operation | Time | Why |
|---|---|---|
| Sort chronologically | O(n log n) | comparison sort |
| Coordinate compression | O(n log n) | sort distinct values; dict lookups are O(1) |
| Build Wavelet Tree | O(n log σ) | each of ~log σ levels touches ≤ n items once |
| Date → index range | O(log n) | two binary searches |
| k-th smallest / largest | O(log σ) | one step (O(1)) per level |
| Count ≥ threshold | O(log σ) | one `bisect` on distinct values + one descent |
| Median | O(log σ) | one or two k-th smallest queries |
| Distribution (4 buckets) | O(log σ) | a constant number of threshold counts |
| Compare two periods | O(log n + log σ) | constant number of the above, ×2 |

**Space:** O(n log σ) integers (each level stores ≤ n prefix entries) + O(n) for the records.
(A production wavelet tree uses bit-vectors with rank structures to get n log σ *bits*; in Python we use plain integer lists for clarity, which costs more real memory.)

## 6. Why a Wavelet Tree and not "sort the selected slice"?

For a slice of length m, sorting costs **O(m log m)** *per query* (or O(m) with quickselect, still per query). A year of USGS data can be thousands of events; a decade, hundreds of thousands. A Wavelet Tree pays **O(n log σ) once**, then each query on any slice costs only **O(log σ)** (about 7–10 steps for magnitudes with 0.1 resolution) *regardless of how big the slice is*. It handles all queries (k-th, rank/count, median) with a single structure, and it needs no per-range preprocessing – unlike precomputed sorted lists or prefix-sum tables, which cannot answer arbitrary range order-statistics. Typical use: an analyst or dashboard asks hundreds of different date-range questions on a fixed catalogue.

## 7. Test cases (all pass: `Ran 18 tests ... OK`)

| Test | What it proves |
|---|---|
| Compression order & duplicates | ranks preserve order, equal values share a rank |
| 150 random arrays × 20 random ranges × every k | k-th smallest/largest equal `sorted()` answers (heavy duplicates) |
| Thresholds incl. −1, 4.95, 100 | ≥-count equals a simple loop count |
| All-equal array, single element | degenerate trees |
| Sample file: validation report | 3 bad magnitudes, 1 bad time, 1 duplicate id dropped |
| Chronological order, one-day range inclusive | sorting + binary-search boundaries |
| Boundary magnitudes 3.0 / 5.0 / 7.0 | each record in exactly one distribution bucket |
| Median odd & even | average of two middle values |
| Empty range, k=0, k>count, start>end | proper errors / `None` instead of crashes |
| User date parsing | `2020`, invalid month/day, wrong format |
| Compare two periods | counts, maximum, threshold counts |
| 20,000 random records × 300 ranges | still correct at scale |

## 8. Limitations

- **Static:** the tree is built once; adding new earthquakes means rebuilding (O(n log σ)).
- **Time resolution is by date:** queries use whole days (inclusive); no time-of-day filtering.
- **Only magnitude is indexed:** queries combining region/depth filters need other structures.
- **Python overhead:** lists of ints use far more memory than the theoretical O(n log σ) *bits*; a 1M-row catalogue is feasible but heavy. Build time is Python-slow compared to C++.
- **Float thresholds** compare exactly to the catalogue's reported decimals (USGS uses 1–2 decimals), so `4.95` vs `4.9` behave as the stored values dictate.
- Duplicate events are detected by `id` only; magnitude corrections/revisions by USGS are not tracked.
- Uses one global magnitude type (`mag`); USGS mixes `mw`, `mb`, `ml`, … which are not strictly comparable.

## 9. Likely viva questions

1. **Why must the data be sorted by time first?** So a date range is one contiguous index range.
2. **Why compress coordinates?** Tree levels split on integer value ranges; floats can't be split by halving reliably, and compression bounds the depth by log σ rather than the value range.
3. **What does `left_prefix` store and why is it enough?** Number of items sent left so far; it lets us translate a range into the child's range in O(1) (see §4).
4. **Why does the tree depth not depend on n?** It halves the *value range*, not the array; depth = ⌈log₂ σ⌉.
5. **How are duplicates handled?** Equal values share one rank and always travel the same path; k-th queries count them as separate elements via range lengths.
6. **Why half-open ranges `[left, right)`?** Length = right − left, empty range is `left == right`, and prefix arrays index naturally.
7. **How is k-th largest derived?** `(length − k + 1)`-th smallest.
"# prakash" 
