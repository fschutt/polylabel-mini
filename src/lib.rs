//! This crate provides a Rust implementation of the [Polylabel](https://github.com/mapbox/polylabel) algorithm
//! for finding the optimum position of a polygon label.
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::f64;

// TODO: polygon.contains(point)
// TODO: point.euclidean_distance(polygon.exterior)
// area()
// centroid()
// bounding_rect()

pub const COORD_PRECISION: f64 = 1e-1; // 0.1m

#[derive(Debug, Copy, Clone)]
pub struct Point {
    pub x: f64,
    pub y: f64,
}

impl PartialEq for Point {
    fn eq(&self, other: &Point) -> bool {
        (self.x - other.x).abs() < COORD_PRECISION &&
        (self.y - other.y).abs() < COORD_PRECISION
    }
}

impl Point {
    fn euclidean_distance_line_string(&self, line_string: &LineString) -> f64 {
        // No need to continue if the point is on the LineString, or it's empty
        if line_string_contains_point(line_string, *self) || line_string.points.is_empty() {
            return 0.0;
        }

        line_string.points.windows(2)
            .map(|line| line_segment_distance(*self, line[0], line[1]))
            .fold(f64::MAX, |accum, val| accum.min(val))
    }
}

fn line_segment_distance(point: Point, start: Point, end: Point) -> f64
{
    if start == end {
        return euclidean_length_2(point, start);
    }
    let dx = dx(start, end);
    let dy = dy(start, end);
    let r = ((point.x - start.x) * dx + (point.y - start.y) * dy) / (dx.powi(2) + dy.powi(2));

    if r <= 0.0 {
        euclidean_length_2(point, start)
    } else if r >= 1.0 {
        euclidean_length_2(point, end)
    } else {
        let s = ((start.y - point.y) * dx - (start.x - point.x) * dy) / (dx * dx + dy * dy);
        s.abs() * dx.hypot(dy)
    }
}

fn line_string_contains_point(line_string: &LineString, point: Point) -> bool {
    // LineString without points
    if line_string.points.is_empty() {
        return false;
    }

    // Check if point is a vertex
    if line_string.points.iter().any(|p| point_contains_point(*p, point)) {
        return true;
    }

    for line in line_string.points.windows(2) {
        let start = line[0];
        let end = line[1];
        if {
            (start.y == end.y)
         && (start.y == point.y)
         && (point.x > start.x.min(end.x))
         && (point.x < start.x.max(end.x))
        } || {
            (start.x == end.x)
         && (start.x == point.x)
         && (point.y > start.y.min(end.y))
         && (point.y < start.y.max(end.y))
        } {
            return true;
        }
    }

    false
}

fn point_contains_point(a: Point, b: Point) -> bool {
    euclidean_length_2(a, b) < COORD_PRECISION
}

#[derive(Debug, Clone)]
pub struct LineString {
    pub points: Vec<Point>,
}

fn centroid_points(start: Point, end: Point) -> Point {
    let two = 2.0;
    let x = start.x + dx(start, end) / two;
    let y = start.y + dy(start, end) / two;
    Point { x, y }
}

impl LineString {
    pub fn area(&self) -> f64 {
        get_linestring_area(&self)
    }

    pub fn centroid(&self) -> Option<Point> {
        // The Centroid of a LineString is the mean of the middle of the segment
        // weighted by the length of the segments.

        if self.points.is_empty() {
            return None;
        }

        if self.points.len() == 1 {
            Some(self.points[0])
        } else {
            let (sum_x, sum_y, total_length) =
                self.points.windows(2)
                    .fold((0.0_f64, 0.0_f64, 0.0_f64), |accum, line| {
                        let segment_len = euclidean_length_2(line[0], line[1]);
                        let line_center = centroid_points(line[0], line[1]);
                        (
                            accum.0 + segment_len * line_center.x,
                            accum.1 + segment_len * line_center.y,
                            accum.2 + segment_len,
                        )
                    });
            Some(Point { x: sum_x / total_length, y: sum_y / total_length })
        }
    }
}

fn dx(a: Point, b: Point) -> f64 {
    b.x - a.x
}

fn dy(a: Point, b: Point) -> f64 {
    b.y - a.y
}

fn euclidean_length_2(a: Point, b: Point) -> f64 {
    dx(a, b).hypot(dy(a, b))
}

#[derive(Debug, Clone)]
pub struct Polygon {
    pub exterior: LineString,
    pub interiors: Vec<LineString>,
}

#[derive(PartialEq, Eq)]
enum PointPosition {
    Inside,
    Outside,
    OnBoundary,
}

impl Polygon {
    fn contains(&self, p: &Point) -> bool {
        match get_position(*p, &self.exterior) {
            PointPosition::OnBoundary | PointPosition::Outside => false,
            _ => self
                .interiors
                .iter()
                .all(|ls| get_position(*p, ls) == PointPosition::Outside),
        }
    }

    fn area(&self) -> f64 {
        let area_exterior = get_linestring_area(&self.exterior);
        let area_interior: f64 = self.interiors.iter().map(|line| get_linestring_area(line)).sum();
        area_exterior - area_interior
    }

    // Calculate the centroid of a Polygon.
    // We distinguish between a simple polygon, which has no interior rings (holes),
    // and a complex polygon, which has one or more interior rings.
    // A complex polygon's centroid is the weighted average of its
    // exterior shell centroid and the centroids of the interior ring(s).
    // Both the shell and the ring(s) are considered simple polygons for the purposes of
    // this calculation.
    //
    // See here for a formula: http://math.stackexchange.com/a/623849
    // See here for detail on alternative methods: https://fotino.me/calculating-centroids/
    fn centroid(&self) -> Option<Point> {

        if self.exterior.points.is_empty() {
            return None;
        } else if self.exterior.points.len() == 1 {
            return Some(self.exterior.points[0]);
        }

        let external_centroid = simple_polygon_centroid(&self.exterior)?;

        if !self.interiors.is_empty() {
            let external_area = simple_polygon_area(&self.exterior).abs();
            // accumulate interior Polygons
            let (totals_x, totals_y, internal_area) = self
                .interiors
                .iter()
                .filter_map(|ring| {
                    let area = simple_polygon_area(ring).abs();
                    let centroid = simple_polygon_centroid(ring)?;
                    Some((centroid.x * area, centroid.y * area, area))
                })
                .fold((0.0_f64, 0.0_f64, 0.0_f64), |accum, val| {
                    (accum.0 + val.0, accum.1 + val.1, accum.2 + val.2)
                });

            return Some(Point {
                x: ((external_centroid.x * external_area) - totals_x) / (external_area - internal_area),
                y: ((external_centroid.y * external_area) - totals_y) / (external_area - internal_area),
            });
        }

        Some(external_centroid)
    }

    fn bounding_rect(&self) -> Option<Rect> {
        get_bounding_rect(&self.exterior)
    }
}


struct Rect {
    pub min: Point,
    pub max: Point,
}

fn get_bounding_rect(line_string: &LineString) -> Option<Rect> {
    if line_string.points.is_empty() {
        return None;
    }
    let first_point = line_string.points[0];

    let mut min_x = first_point.x;
    let mut max_x = first_point.x;
    let mut min_y = first_point.y;
    let mut max_y = first_point.y;

    for point in line_string.points.iter() {
        min_x = min_x.min(point.x);
        max_x = max_x.max(point.x);
        min_y = min_y.min(point.y);
        max_y = max_y.max(point.y);
    }

    Some(Rect {
        min: Point { x: min_x, y: min_y },
        max: Point { x: max_x, y: max_y },
    })
}

fn simple_polygon_area(line_string: &LineString) -> f64 {
    if line_string.points.is_empty() || line_string.points.len() == 1 {
        return 0.0;
    }

    line_string.points.windows(2).map(|line| determinant(line[0], line[1])).sum::<f64>() / 2.0_f64
}

fn get_linestring_area(line_string: &LineString) -> f64 {
    line_string.points.windows(2).map(|line| {
        determinant(line[0], line[1])
    }).sum()
}

fn determinant(start: Point, end: Point) -> f64 {
    start.x * end.y - start.y * end.x
}

fn simple_polygon_centroid(line_string: &LineString) -> Option<Point> {
    let area = get_linestring_area(line_string);
    if area == 0.0 {
        // if the polygon is flat (area = 0), it is considered as a linestring
        return line_string.centroid();
    }
    let (sum_x, sum_y) = line_string.points.windows(2)
        .fold((0.0_f64, 0.0_f64), |accum, line| {
            let start = line[0];
            let end = line[1];
            let tmp = determinant(start, end);
            (
                accum.0 + ((end.x + start.x) * tmp),
                accum.1 + ((end.y + start.y) * tmp),
            )
        });

    let six = 6.0;
    Some(Point { x: sum_x / (six * area), y: sum_y / (six * area) })
}

/// Calculate the position of `Point` p relative to a linestring
fn get_position(p: Point, linestring: &LineString) -> PointPosition {

    // See:
    // http://geospatialpython.com/search?updated-min=2011-01-01T00:00:00-06:00&updated-max=2012-01-01T00:00:00-06:00&max-results=19

    // LineString without points
    if linestring.points.is_empty() {
        return PointPosition::Outside;
    }
    // Point is on linestring
    if line_string_contains_point(&linestring, p) {
        return PointPosition::OnBoundary;
    }

    let mut xints = 0.0_f64;
    let mut crossings = 0;

    for line in linestring.points.windows(2) {
        let start = line[0];
        let end = line[1];
        if p.y > start.y.min(end.y)
            && p.y <= start.y.max(end.y)
            && p.x <= start.x.max(end.x)
        {
            if start.y != end.y {
                xints = (p.y - start.y) * (end.x - start.x)
                    / (end.y - start.y)
                    + start.x;
            }
            if (start.x == end.x) || (p.x <= xints) {
                crossings += 1;
            }
        }
    }

    if crossings % 2 == 1 {
        PointPosition::Inside
    } else {
        PointPosition::Outside
    }
}


/// Represention of a Quadtree node's cells. A node contains four Qcells.
#[derive(Debug)]
struct Qcell {
    // The cell's centroid
    centroid: Point,
    // Half of the parent node's extent
    extent: f64,
    // Distance from centroid to polygon
    distance: f64,
    // Maximum distance to polygon within a cell
    max_distance: f64,
}

impl Qcell {
    fn new(x: f64, y: f64, h: f64, distance: f64, max_distance: f64) -> Qcell {
        Qcell {
            centroid: Point { x, y },
            extent: h,
            distance,
            max_distance,
        }
    }
}

impl PartialOrd for Qcell {
    fn partial_cmp(&self, other: &Qcell) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Qcell {
    fn cmp(&self, other: &Qcell) -> std::cmp::Ordering {
        self.max_distance.partial_cmp(&other.max_distance).unwrap()
    }
}

impl PartialEq for Qcell {
    fn eq(&self, other: &Qcell) -> bool {
        (self.max_distance - other.max_distance).abs() < COORD_PRECISION
    }
}

impl Eq for Qcell { }

/// Signed distance from a Qcell's centroid to a Polygon's outline
/// Returned value is negative if the point is outside the polygon's exterior ring
fn signed_distance(x: f64, y: f64, polygon: &Polygon) -> f64 {
    let point = Point { x, y };
    let inside = polygon.contains(&point);
    // Use LineString distance, because Polygon distance returns 0.0 for inside
    let distance = point.euclidean_distance_line_string(&polygon.exterior);
    if inside {
        distance
    } else {
        -distance
    }
}

/// Add a new Quadtree node made up of four `Qcell`s to the binary heap
fn add_quad(
    mpq: &mut BinaryHeap<Qcell>,
    cell: &Qcell,
    new_height: f64,
    polygon: &Polygon,
) {
    let two = 2.0_f64;
    let centroid_x = cell.centroid.x;
    let centroid_y = cell.centroid.y;
    for combo in &[
        (centroid_x - new_height, centroid_y - new_height),
        (centroid_x + new_height, centroid_y - new_height),
        (centroid_x - new_height, centroid_y + new_height),
        (centroid_x + new_height, centroid_y + new_height),
    ] {
        let mut new_dist = signed_distance(combo.0, combo.1, polygon);
        mpq.push(Qcell::new(
            combo.0,
            combo.1,
            new_height,
            new_dist,
            new_dist + new_height * two.sqrt(),
        ));
    }
}

/// Calculate a Polygon's ideal label position by calculating its ✨pole of inaccessibility✨
///
/// The calculation uses an [iterative grid-based algorithm](https://github.com/mapbox/polylabel#how-the-algorithm-works).
///
/// # Examples
///
/// ```
/// use polylabel::polylabel;
/// extern crate geo;
/// use self::geo::{Point, LineString, Polygon};
///
/// // An approximate `L` shape
/// let coords = vec![
///    (0.0, 0.0),
///    (4.0, 0.0),
///    (4.0, 1.0),
///    (1.0, 1.0),
///    (1.0, 4.0),
///    (0.0, 4.0),
///    (0.0, 0.0)];
///
/// let poly = Polygon::new(coords.into(), vec![]);
///
/// // Its centroid lies outside the polygon
/// assert_eq!(poly.centroid(), Point::new(1.3571428571428572, 1.3571428571428572));
///
/// let label_position = polylabel(&poly, &1.0);
/// // Optimum label position is inside the polygon
/// assert_eq!(label_position, Point::new(0.5625, 0.5625));
/// ```
///
pub fn polylabel(polygon: &Polygon, tolerance: f64) -> Point {

    // special case for degenerate polygons
    if polygon.area() < 0.0 {
        return Point { x: 0.0, y: 0.0 };
    }

    let two = 2.0_f64;

    // Initial best cell values
    let centroid = polygon.centroid().unwrap();
    let bbox = polygon.bounding_rect().unwrap();
    let width = bbox.max.x - bbox.min.x;
    let height = bbox.max.y - bbox.min.y;
    let cell_size = width.min(height);

    // Special case for degenerate polygons
    if cell_size == 0.0 {
        return Point { x: bbox.min.x, y: bbox.min.y };
    }

    let mut h = cell_size / two;
    let distance = signed_distance(centroid.x, centroid.y, polygon);
    let max_distance = distance + 0.0 * two.sqrt();

    let mut best_cell = Qcell::new(
        centroid.x,
        centroid.y,
        0.0,
        distance,
        max_distance,
    );

    // special case for rectangular polygons
    let bbox_cell_dist = signed_distance(
        bbox.min.x + width / two,
        bbox.min.y + height / two,
        polygon,
    );
    let bbox_cell = Qcell {
        centroid: Point { x: bbox.min.x + width / two, y: bbox.min.y + height / two },
        extent: 0.0,
        distance: bbox_cell_dist,
        max_distance: bbox_cell_dist + 0.0 * two.sqrt(),
    };

    if bbox_cell.distance > best_cell.distance {
        best_cell = bbox_cell;
    }

    // Priority queue
    let mut cell_queue: BinaryHeap<Qcell> = BinaryHeap::new();
    // Build an initial quadtree node, which covers the Polygon
    let mut x = bbox.min.x;
    let mut y;

    while x < bbox.max.x {
        y = bbox.min.y;
        while y < bbox.max.y {
            let latest_dist = signed_distance(x + h, y + h, polygon);
            cell_queue.push(Qcell {
                centroid: Point { x: x + h, y: y + h },
                extent: h,
                distance: latest_dist,
                max_distance: latest_dist + h * two.sqrt(),
            });
            y = y + cell_size;
        }
        x = x + cell_size;
    }


    // Now try to find better solutions
    while !cell_queue.is_empty() {
        let cell = cell_queue.pop().unwrap();
        // Update the best cell if we find a cell with greater distance
        if cell.distance > best_cell.distance {
            best_cell.centroid = Point { x: cell.centroid.x, y: cell.centroid.y };
            best_cell.extent = cell.extent;
            best_cell.distance = cell.distance;
            best_cell.max_distance = cell.max_distance;
        }
        // Bail out of this iteration if we can't find a better solution
        if cell.max_distance - best_cell.distance <= tolerance {
            continue;
        }
        // Otherwise, add a new quadtree node and start again
        h = cell.extent / two;
        add_quad(&mut cell_queue, &cell, h, polygon);
    }

    // We've exhausted the queue, so return the best solution we've found
    best_cell.centroid
}