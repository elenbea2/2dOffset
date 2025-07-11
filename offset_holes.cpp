#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h> 
#include <CGAL/intersections.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <cmath>
#include <CGAL/Polygon_2_algorithms.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Vector_2 Vector;
typedef K::Segment_2 Segment;
typedef K::Line_2 Line;

struct JumpLink {
    size_t index_a;
    size_t index_b;
};


double signed_area(const std::vector<Point>& poly) {
    double area = 0.0;
    int n = poly.size();
    for (int i = 0; i < n; ++i) {
        const Point& p1 = poly[i];
        const Point& p2 = poly[(i + 1) % n];
        area += (p1.x() * p2.y()) - (p2.x() * p1.y());
    }
    return area / 2.0;
}

std::vector<std::vector<Point>> remove_clockwise_polygons(
    const std::vector<std::vector<Point>>& polygons,
    const std::set<std::pair<double, double>>& intersection_points)
{
    std::vector<std::vector<Point>> result;
    for (size_t i = 0; i < polygons.size(); ++i) {
        double area = signed_area(polygons[i]);
        std::cout << "→ Polígono " << i << " | Área assinada: " << area << "\n";
        if (area > 0) {
            result.push_back(polygons[i]);
        } else {
            bool has_intersection = false;
            for (const auto& p : polygons[i]) {
                if (intersection_points.count({p.x(), p.y()})) {
                    has_intersection = true;
                    std::cout << "    Buraco removido: ponto (" << p.x() << ", " << p.y() << ") está em interseção\n";
                    break;
                }
            }
            if (!has_intersection) {
                std::cout << "    Buraco mantido (sem interseção)\n";
                result.push_back(polygons[i]);
            }
        }
    }
    return result;
}

std::vector<std::pair<Point, std::pair<std::pair<int, int>, std::pair<int, int>>>> 
find_all_self_intersections(const std::vector<std::vector<Point>>& polygons) {
    std::vector<std::pair<Point, std::pair<std::pair<int, int>, std::pair<int, int>>>> intersections;
    int poly_count = polygons.size();

    for (int pi = 0; pi < poly_count; ++pi) {
        const auto& polyA = polygons[pi];
        int nA = polyA.size();
        for (int i = 0; i < nA; ++i) {
            Segment s1(polyA[i], polyA[(i + 1) % nA]);

            for (int pj = pi; pj < poly_count; ++pj) {
                const auto& polyB = polygons[pj];
                int nB = polyB.size();
                int j_start = (pi == pj) ? i + 2 : 0;

                for (int j = j_start; j < nB; ++j) {
                    if (pi == pj && i == 0 && j == nB - 1) continue;

                    Segment s2(polyB[j], polyB[(j + 1) % nB]);
                    auto res = CGAL::intersection(s1, s2);
                    if (res) {
                        if (const Point* p = boost::get<Point>(&*res)) {
                            intersections.push_back({*p, {{pi, i}, {pj, j}}});
                        }
                    }
                }
            }
        }
    }
    return intersections;
}

double param_along_segment(const Point& a, const Point& b, const Point& p) {
    Vector ab = b - a;
    Vector ap = p - a;
    return (ab * ap) / ab.squared_length(); 
}


void insert_intersections_and_links_multi(
    std::vector<std::vector<Point>>& polygons,
    const std::vector<std::pair<Point, std::pair<std::pair<int, int>, std::pair<int, int>>>>& intersections,
    std::vector<JumpLink>& links,
    std::vector<std::vector<size_t>>& poly_map 
) {

    std::vector<std::vector<Point>> new_polys(polygons.size());
    std::vector<std::map<int, std::vector<std::pair<Point, int>>>> inserts(polygons.size());
    std::vector<size_t> copy_indices_a(intersections.size());
    std::vector<size_t> copy_indices_b(intersections.size());

    for (size_t k = 0; k < intersections.size(); ++k) {
        const auto& inter = intersections[k];
        int pi_a = inter.second.first.first;
        int seg_a = inter.second.first.second;
        int pi_b = inter.second.second.first;
        int seg_b = inter.second.second.second;
        inserts[pi_a][seg_a].push_back({inter.first, static_cast<int>(k * 2)});
        inserts[pi_b][seg_b].push_back({inter.first, static_cast<int>(k * 2 + 1)});
    }

    for (size_t p = 0; p < polygons.size(); ++p) {
        for (auto& [seg_idx, points] : inserts[p]) {
            const Point& a = polygons[p][seg_idx];
            const Point& b = polygons[p][(seg_idx + 1) % polygons[p].size()];
            std::sort(points.begin(), points.end(),
                [&](const std::pair<Point, int>& lhs, const std::pair<Point, int>& rhs) {
                    return param_along_segment(a, b, lhs.first) < param_along_segment(a, b, rhs.first);
                });
        }
    }
    
    size_t global_idx = 0;
    poly_map.clear();
    poly_map.resize(polygons.size());
    for (size_t p = 0; p < polygons.size(); ++p) {
        auto& poly = polygons[p];
        auto& new_poly = new_polys[p];
        auto& map_vec = poly_map[p];
        int count = 0;
        std::cout << "Polígono " << p << ":\n";
        for (size_t i = 0; i < poly.size(); ++i) {
            new_poly.push_back(poly[i]);
            map_vec.push_back(global_idx++);
            std::cout << "  Ponto original: (" << poly[i].x() << ", " << poly[i].y() << ")\n";
            count++;
            if (inserts[p].count(i)) {
                for (auto& ins : inserts[p][i]) {
                    new_poly.push_back(ins.first);
                    map_vec.push_back(global_idx++);
                    std::cout << "  Interseção inserida: (" << ins.first.x() << ", " << ins.first.y() << ")\n";
                    int k = ins.second / 2;
                    if (ins.second % 2 == 0)
                        copy_indices_a[k] = map_vec.back();
                    else
                        copy_indices_b[k] = map_vec.back();
                    count++;
                }
            }
        }
    }

    std::set<size_t> used_indices;
    for (size_t k = 0; k < intersections.size(); ++k) {
        size_t idx_a = copy_indices_a[k];
        size_t idx_b = copy_indices_b[k];
        if (used_indices.find(idx_a) == used_indices.end() &&
            used_indices.find(idx_b) == used_indices.end()) {
            links.push_back({idx_a, idx_b});
            used_indices.insert(idx_a);
            used_indices.insert(idx_b);
        } else {
            std::cerr << "⚠️ Problema ao criar link multi " << k << ": índices já usados.\n";
        }
    }

    
    std::cout << "\nLINKS\n";
    for (size_t k = 0; k < links.size(); ++k) {
        size_t idx_a = links[k].index_a;
        size_t idx_b = links[k].index_b;
        std::cout << "Link " << k << ": idx_a = " << idx_a << " <--> idx_b = " << idx_b << "\n";
    }

    polygons = new_polys;
}

std::vector<std::vector<Point>> split_polygons_multi(
    const std::vector<std::vector<Point>>& polygons,
    const std::vector<JumpLink>& links)
{
    std::map<size_t, std::pair<size_t, size_t>> global_to_local;
    std::map<std::pair<size_t, size_t>, size_t> local_to_global;

    size_t global_idx = 0;
    for (size_t pi = 0; pi < polygons.size(); ++pi) {
        for (size_t li = 0; li < polygons[pi].size(); ++li) {
            global_to_local[global_idx] = {pi, li};
            local_to_global[{pi, li}] = global_idx;
            ++global_idx;
        }
    }

    std::map<size_t, size_t> jump_map;
    for (const auto& link : links) {
        jump_map[link.index_a] = link.index_b;
        jump_map[link.index_b] = link.index_a;
    }

    std::set<size_t> visited;
    std::vector<std::vector<Point>> fragments;

    for (const auto& [start_global, _] : global_to_local) {
        if (visited.count(start_global)) continue;

        std::vector<Point> frag;
        size_t current = start_global;
        bool first = true;
        auto [current_pi, current_li] = global_to_local[current];

        while (first || (!visited.count(current))) {
            first = false;

            frag.push_back(polygons[current_pi][current_li]);
            visited.insert(current);

            if (jump_map.count(current)) {
                size_t target = jump_map[current];
                if (!global_to_local.count(target)) break;

                auto [pi, li] = global_to_local[target]; 
                size_t next_li = (li + 1) % polygons[pi].size();
                size_t next = local_to_global[{pi, next_li}];

                if (visited.count(next)) break;

                current = next;
                std::tie(current_pi, current_li) = global_to_local[current];
                continue;

            }

            size_t next_li = (current_li + 1) % polygons[current_pi].size();
            size_t next_global = local_to_global[{current_pi, next_li}];

            if (visited.count(next_global)) break;

            current_li = next_li;
            current = next_global;
        }

        if (frag.size() >= 1)
            fragments.push_back(frag);
    }

    
    std::cout << "\nFRAGMENTOS (índices globais)\n";
    for (size_t i = 0; i < fragments.size(); ++i) {
        std::cout << "Fragmento " << i << " (" << fragments[i].size() << " pontos): ";
        for (const auto& p : fragments[i]) {
            
            for (const auto& [idx, local] : global_to_local) {
                if (polygons[local.first][local.second] == p) {
                    std::cout << idx << " ";
                    break;
                }
            }
        }
        std::cout << "\n";
    }

    return fragments;
}



Line offset_edge(const Point& p1, const Point& p2, double d) {
    Vector v = p2 - p1;
    Vector normal(v.y(), -v.x());
    normal = normal / std::sqrt(normal.squared_length());
    Point q1 = p1 + normal * d;
    Point q2 = p2 + normal * d;
    return Line(q1, q2);
}

std::vector<Line> offset_polygon(const std::vector<Point>& poly, const std::vector<double>& distances) {
    std::vector<Line> offset_lines;
    int n = poly.size();
    for (int i = 0; i < n; ++i) {
        const Point& p1 = poly[i];
        const Point& p2 = poly[(i+1)%n];
        offset_lines.push_back(offset_edge(p1, p2, distances[i]));
    }
    return offset_lines;
}



std::vector<Point> compute_edge_connections(const std::vector<Line>& lines) {
    std::vector<Point> intersections;
    int n = lines.size();
    for (int i = 0; i < n; ++i) {
        auto res = CGAL::intersection(lines[i], lines[(i+1)%n]);
        if (res) {
            if (const Point* p = boost::get<Point>(&*res)) {
                intersections.push_back(*p);
            } else {
                std::cerr << "Aviso: interseção não é um ponto entre linhas " << i << " e " << (i+1)%n << std::endl;
            }
        } else {
            std::cerr << "Aviso: sem interseção visível entre linhas " << i << " e " << (i+1)%n << std::endl;
        }
    }
    return intersections;
}

void export_svg(const std::vector<std::vector<Point>>& outers, 
                const std::vector<std::vector<Point>>& sub_polygons, 
                const std::string& filename,
                const std::vector<std::vector<Point>>& holes = {})
{
    std::ofstream out(filename);
    out << "<svg xmlns='http://www.w3.org/2000/svg' width='600' height='600' viewBox='0 -60 100 140'>\n";

    for (const auto& original : outers) {
        out << "<polygon points='";
        for (auto& p : original)
            out << p.x() << "," << -p.y() << " ";
        out << "' fill='none' stroke='blue' stroke-width='2'/>\n";
    }

    for (const auto& hole : holes) {
        out << "<polygon points='";
        for (auto& p : hole)
            out << p.x() << "," << -p.y() << " ";
        out << "' fill='none' stroke='blue' stroke-width='2' stroke-dasharray='4,2'/>\n";
    }

    for (const auto& deslocado : sub_polygons) {
        out << "<polygon points='";
        for (auto& p : deslocado)
            out << p.x() << "," << -p.y() << " ";
        out << "' fill='none' stroke='red' stroke-width='1'/>\n";

        
        for (auto& p : deslocado)
            out << "<circle cx='" << p.x() << "' cy='" << -p.y() << "' r='0.5' fill='red'/>\n";
    }

    out << "</svg>\n";
    out.close();
}

int main() {
    
    std::vector<Point> hole1 = {
        Point(5, 35),
        Point(15, 35),
        Point(15, 25),
        Point(5, 25)
    };

    std::vector<Point> hole2 = {
        Point(25, 15),
        Point(35, 15),
        Point(35, 5),
        Point(25, 5) 
    };

    std::vector<std::vector<Point>> holes = {hole1, hole2};

    std::vector<std::vector<double>> distances_holes = {
        {-5, -5, -5, -5},
        {-4, -4, -4, -4}
    };

    
    std::vector<std::vector<Point>> offset_holes_points;
    for (size_t i = 0; i < holes.size(); ++i) {
        const auto& h = holes[i];
        const auto& dist = (i < distances_holes.size()) ? distances_holes[i] : std::vector<double>(h.size(), -2);

        
        auto offset_hole_lines = offset_polygon(h, dist);
        auto offset_hole_points = compute_edge_connections(offset_hole_lines);

        offset_holes_points.push_back(offset_hole_points);
    }

    
    std::vector<Point> outer1 = {
        Point(0, 0),
        Point(40, 0),
        Point(40, 40),
        Point(0, 40)
    };    
    
    std::vector<Point> outer2 = {
        Point(60, 0),
        Point(100, 0),
        Point(100, 40),
        Point(60, 40)
    };
    
    /*std::vector<Point> outer1 = {
        Point(0, 10),
        Point(5, -10),
        Point(-5, -5),
        Point(-5, -20),
        Point(5, -15),
        Point(0, -40),
        Point(20, -40),
        Point(15, -15),
        Point(35, -20),
        Point(35, -5),
        Point(15, -10),
        Point(20, 10)
    };*/

    /*std::vector<Point> outer1 = {
        Point(60, 40),
        Point(20, 40),
        Point(30, 20),
        Point(10, -10),
        Point(40, -30),
        Point(50, 10),
        Point(70, -20),
        Point(80, 20),
        Point(50, 20)
    };*/
    

    std::vector<std::vector<Point>> outers = {outer1, outer2};

    std::vector<std::vector<double>> distances_outers = {
        //{-4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0}
        //{-10.5 -10.5, -10.5, -10.5, -10.5, -10.5, -10.5, -10.5, -10.5}
        {-2, -2, -2, -2},
        {-2, -2, -2, -2}
    };

    
    std::vector<std::vector<Point>> offset_outers_points;
    for (size_t i = 0; i < outers.size(); ++i) {
        const auto& o = outers[i];
        const auto& dist = (i < distances_outers.size()) ? distances_outers[i] : std::vector<double>(o.size(), -2);
        auto offset_lines = offset_polygon(o, dist);
        auto intersections = compute_edge_connections(offset_lines);

        
        std::cout << "Contorno externo " << i << ":\n";
        for (size_t j = 0; j < o.size(); ++j) {
            std::cout << "  Original: (" << o[j].x() << ", " << o[j].y() << ")"
                      << "  ->  Offsetado: (" << intersections[j].x() << ", " << intersections[j].y() << ")\n";
        }

        offset_outers_points.push_back(intersections);
    }
    
    std::vector<std::vector<Point>> all_offset_polys;
    for (const auto& o : offset_outers_points)
        all_offset_polys.push_back(o);
    for (const auto& h : offset_holes_points)
        all_offset_polys.push_back(h);

        
    auto all_inters = find_all_self_intersections(all_offset_polys);

    
    std::set<std::pair<double, double>> intersection_points;
    for (const auto& inter : all_inters) {
        intersection_points.insert({inter.first.x(), inter.first.y()});
    }

    std::vector<JumpLink> links;
    std::vector<std::vector<size_t>> poly_map;
    insert_intersections_and_links_multi(all_offset_polys, all_inters, links, poly_map);

    auto sub_polygons = split_polygons_multi(all_offset_polys, links);

    
    std::set<size_t> link_indices;
    for (const auto& link : links) {
        link_indices.insert(link.index_a);
        link_indices.insert(link.index_b);
    }

    
    std::vector<bool> polygon_has_intersection(sub_polygons.size(), false);
    size_t global_idx = 0;
    for (size_t i = 0; i < sub_polygons.size(); ++i) {
        bool found = false;
            for (size_t j = 0; j < sub_polygons[i].size(); ++j) {
                const auto& p = sub_polygons[i][j];
                if (intersection_points.count({p.x(), p.y()})) {
                    found = true;
                    std::cout << "Polígono " << i << " contém interseção em ponto " << j
                            << " (" << p.x() << ", " << p.y() << ")\n";
                    break;
                }
            }

        polygon_has_intersection[i] = found;
    }

    

    auto filtered_polygons = remove_clockwise_polygons(sub_polygons, intersection_points);

    
    export_svg(outers, filtered_polygons, "output.svg", holes);

    std::cout << "Arquivo 'output.svg' gerado.\n";

    return 0;
}
