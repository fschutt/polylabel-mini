# polylabel-mini

This is a minimal implementation of the polylabel algorithm, for `f64` polygons
(since usually, that's what you have to work with in GIS applications anyways, single-precision is too unprecise).
Original code can be found at: https://github.com/urschrei/polylabel-rs

The reason for this fork was that `polylabel-rs` requires 18 (!) dependencies for
an algorithm that is merely 50 lines long. This leads to a huge bloat in compile times.

Yes, this is technically a copy-paste implementation, but for my purposes, I only
use the `polylabel()` function and nothing else. It simply doesn't make sense to
import 18 extremely generics-heavy dependencies for 50 lines of code.

This crate doesn't have any dependencies and has special optimizations applied, for example
calculating the bounding rectangle doesn't require cloning the points.

The following functions were copied from the `geo` crate:

```
- Polygon::contains(&Point)
- Point::euclidean_distance(&Polygon)
- Polygon::area()
- Polygon::centroid()
- Polygon::bounding_rect()
```

As you can see, not anything that would justify 18 dependencies. The code was forked at
`master#fc942fab47ca29cbc54597448b7af48d7cad109c` and is updated from time to time.

The license of both `geo` and `polylabel` is MIT, therefore this repository is also MIT
licensed. I do not want to violate any licensing restrictions, I just want to reduce
dependencies. All code is licensed under the MIT license from Stefan Hügel and Corey Farwell.
