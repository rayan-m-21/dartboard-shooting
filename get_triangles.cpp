/**
 * get_triangles.cpp
 *
 * Implementation of polygon triangulation using the Ear Clipping algorithm.
 *
 * This program reads a list of 2D points describing a simple polygon
 * (assumed to be given in counter-clockwise order) and outputs a set of
 * non-overlapping triangles that fully cover the polygon.
 *
 * The algorithm works by iteratively identifying and removing "ears"
 * (triangles formed by a convex vertex and its neighbors that contain no
 * other vertices inside). Each clipped ear is recorded until the polygon
 * is reduced to a single triangle.
 *
 * ----------------------------------------------------------------------
 * Command Line Usage:
 *
 *   ./get_triangles <polygon_points_file> <outfile>
 *
 * Arguments:
 *   <polygon_points_file> : A text file containing the polygon vertices,
 *                           listed one per line in the format:
 *                              x y
 *                           Vertices must be provided in counter-clockwise
 *                           order. The polygon must be simple (non-self
 *                           intersecting).
 *
 *   <outfile>             : Path to the file where the resulting triangles
 *                           will be written. Each line represents one
 *                           triangle, written as three vertices separated
 *                           by "|" characters:
 *                              x1 y1|x2 y2|x3 y3
 *
 * Example:
 *   Input file (square.txt):
 *       0 0
 *       1 0
 *       1 1
 *       0 1
 *
 *   Command:
 *       ./get_triangles square.txt triangles.txt
 *
 *   Output file (triangles.txt):
 *       0 0|1 0|1 1
 *       0 0|1 1|0 1
 *
 * ----------------------------------------------------------------------
 * Notes:
 *  - Debugging mode can be enabled by setting `dbg = true;` to print the
 *    internal state during triangulation.
 *
 */

#include <iostream>
#include <fstream>
using namespace std;

#include <vector>
#include <queue>
#include <set>
#include <cassert>

const double EPS = 1e-9;

bool dbg = false;
string input_file;

#define pdd pair<double, double>

vector<pair<double, double> > points;
vector<int> next_v; // holds the vertex which is currently next on the polygon (anticlock)
vector<int> prev_v;
set<int> active_vertices; // labelled by their index
set<int> active_reflex;
set<int> active_convex;
set<int> active_ears;
set<tuple<int, int, int>> triangles_found; // stored as their indices

void print_state() {
    cout << "___________________________________\n";
    cout << "Active Vertices:\n";
    for (int v: active_vertices) {
        cout << v << " ";
    }
    cout << "\n";

    cout << "reflex:\n";
    for (int i: active_reflex) {
        cout << i << " ";
    }
    cout << '\n';

    cout << "ears:\n";
    for (int i: active_ears) {
        cout << i << " ";
    }
    cout << '\n';

    cout << "convex:\n";
    for (int i: active_convex) {
        cout << i << " ";
    }
    cout << '\n';

    cout << "Found Triangles:\n";
    for (const auto& tri: triangles_found) {
        cout << get<0>(tri) << " "<< get<1>(tri) << " "<< get<2>(tri) << "\n";
    }

    for (int x: active_ears) {
        if (active_convex.find(x) == active_convex.end()) {
            cout << "ERROR: " << x << " is an ear but not a convex?!\n";
            exit(1);
        }
    }

    cout << "DONE PRINTING DATA\n";
    cout << "___________________________________\n";
}

double cross(pdd a, pdd b) {
    return a.first*b.second - a.second*b.first;
}

pdd operator-(pdd a, pdd b) {
  return pdd(a.first - b.first, a.second - b.second);
}

pdd get_prev(int i) {
    if (i >= prev_v.size()) {
        cout << "OUT OF INDEX";
        exit(1);
    }

    return points[prev_v[i]];
}

pdd get_next(int i) {
    if (i >= next_v.size()) {
        cout << "OUT OF INDEX";
        exit(1);
    }

    return points[next_v[i]];
}

bool points_equal(pdd& p1, pdd& p2) {
    double dx = p1.first - p2.first;
    double dy = p1.second - p2.second;
    return (dx * dx + dy * dy) < EPS;
}

double double_area(pdd x1, pdd x2, pdd x3) {
    return abs(cross(x2 - x1, x3 - x1));
}

bool in_triangle(pdd p1, pdd p2, pdd p3, pdd point) {
    double area_triangle = double_area(p1, p2, p3);
    double area1 = double_area(point, p2, p3);
    double area2 = double_area(p1, point, p3);
    double area3 = double_area(p1, p2, point);

    return (abs(area_triangle - (area1 + area2 + area3)) < EPS);
}

bool is_reflex(int i) {
    if (points.size() < 3) exit(1);

    pdd prev = get_prev(i);
    pdd curr = points[i];
    pdd next = get_next(i);

    // if it makes a clockwise turn then it is a reflex point
    return (cross(curr - prev, next - prev) <= EPS);
}

bool is_convex(int i) {
    return !is_reflex(i);
}

set<int> get_reflex_points() {
    set<int> ret = set<int>();

    for (int i = 0; i < points.size(); i++) {
        if (is_reflex(i)) ret.insert(i);
    }

    return ret;
}

set<int> get_convex_points() {
    set<int> ret = set<int>();

    for (int i = 0; i < points.size(); i++) {
        if (!is_reflex(i)) ret.insert(i);
    }

    return ret;
}

bool is_ear(int i) {
    if (is_reflex(i)) return false;

    pdd prev = get_prev(i);
    pdd curr = points[i];
    pdd next = get_next(i);

    for (int j: active_reflex) {
        if (points_equal(points[j], prev) || points_equal(points[j], next)) {
            continue;
        }

        if (in_triangle(prev, curr, next, points[j])) {
            return false;
        }
    }

    return true;
}

set<int> get_ears() {
    // a vertex is an ear iff it is convex and doesnt contain any
    // other vertices in it's triangle
    // Note that it is sufficient to only check that it doesnt contain
    // any of the reflex vertices
    set<int> ret = set<int>();

    for (int i: active_convex) {
        if (is_ear(i)) ret.insert(i);
    }

    return ret;
}

void remove_collinear_points() {
    vector<pdd> new_points;

    for (int i = 0; i < points.size(); i++) {
        pdd prev = get_prev(i);
        pdd curr = points[i];
        pdd next = get_next(i);

        if (abs(cross(curr - prev, next - prev)) > EPS) {
            new_points.push_back(curr);
        }
    }

    points = new_points;
}

/**
 * Will update the status of the given vertex
 */
void update_vertex_status(int i) {
    if (!is_convex(i) && active_convex.find(i) != active_convex.end()) {
        // CASE: No longer a convex but is in convex set
        // So it is a reflex
        active_convex.erase(i);
        active_reflex.insert(i);
    } else if (!is_reflex(i) && active_reflex.find(i) != active_reflex.end()) {
        // CASE: No longer a reflex
        // So it is a convex
        active_convex.insert(i);
        active_reflex.erase(i);
    }

    // cerr << "Checking if " << i << " might be an ear?\n";
    // vertex might be an ear
    if (is_ear(i)) {
        active_ears.insert(i);
    } else {
        active_ears.erase(i);
    }
}

/**
 * This function recursively implements the clipping of ears algorithm
 */
void clip_ears() {
    if (dbg) print_state();

    int num_verts_init = active_vertices.size();

    if (num_verts_init == 3) {
        // can safely add in the vertices remaining
        auto it = active_vertices.begin();
        triangles_found.insert({*it, *next(it), *next(next(it))});
        return;
    }

    if (active_ears.empty()) {
        cout << "ALGORITHM FAILED on " << input_file << "\n";
        exit(1);
    }

    int curr_ind = *active_ears.begin();

    // Remove curr from all sets
    active_vertices.erase(curr_ind);
    active_convex.erase(curr_ind);
    active_ears.erase(curr_ind);

    // Get next and previos indices
    if (dbg) cout << "REMOVING " << curr_ind << "\n";
    int next_ind = next_v[curr_ind];
    int prev_ind = prev_v[curr_ind];

    // update next, prev references
    next_v[prev_ind] = next_ind;
    prev_v[next_ind] = prev_ind;

    // update the status (wether reflex, convex, ear)
    update_vertex_status(next_ind);
    update_vertex_status(prev_ind);

    assert(num_verts_init = active_vertices.size() + 1);

    // add triangle
    triangles_found.insert({curr_ind, prev_ind, next_ind});

    // Clip an ear again to continue the recursion
    clip_ears();
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        cout << "Usage: " << "get_triangles <polygon_points_file> <outfile>\n";
        exit(1);
    }

    input_file = argv[1];

    // Read from the text file
    ifstream infile(argv[1]);

    if (!infile) {
        cerr << "Failed to open file.\n";
        return 1;
    }

    double x, y;
    while (infile >> x >> y) {
        points.push_back(make_pair(x, y));
    }

    // Close the file
    infile.close();

    for (int i = 0; i < points.size(); i++) {
        active_vertices.insert(i);
        next_v.push_back(i + 1);
        prev_v.push_back(i - 1);
    }

    // remove_collinear_points();

    next_v[points.size() - 1] = 0;
    prev_v[0] = points.size() - 1;

    active_convex = get_convex_points();
    active_reflex = get_reflex_points();
    active_ears = get_ears();

    clip_ears();

    if (dbg) print_state();

    ofstream outfile(argv[2]);

    if (!outfile) {
        cerr << "Error opening file: " << argv[2] << "\n";
        return 1;
    }

    for (const auto& tri: triangles_found) {
        int p1= get<0>(tri);
        int p2= get<1>(tri);
        int p3= get<2>(tri);
        outfile << points[p1].first << " " << points[p1].second << "|";
        outfile << points[p2].first << " " << points[p2].second << "|";
        outfile << points[p3].first << " " << points[p3].second << "\n";
    }
}