#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/intersections.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <cmath>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Vector_2 Vector;
typedef K::Segment_2 Segment;
typedef K::Line_2 Line;
typedef CGAL::Polygon_2<K> Polygon;

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

std::vector<std::vector<Point>> remove_clockwise_polygons(const std::vector<std::vector<Point>>& polygons) {
    std::vector<std::vector<Point>> result;
    for (const auto& poly : polygons) {
        double area = signed_area(poly);
        if (area > 0) {
            result.push_back(poly);  // Mantém apenas os anti-horários
        }
    }
    return result;
}

// Função para detectar interseções internas
std::vector<std::pair<Point, std::pair<int, int>>> find_self_intersections(const std::vector<Point>& polygon) {
    std::vector<std::pair<Point, std::pair<int, int>>> intersections;
    int n = polygon.size();
    for (int i = 0; i < n; ++i) {
        Segment s1(polygon[i], polygon[(i + 1) % n]);
        for (int j = i + 2; j < n; ++j) {
            if ((i == 0 && j == n - 1)) continue;
            Segment s2(polygon[j], polygon[(j + 1) % n]);
            auto res = CGAL::intersection(s1, s2);
            if (res) {
                if (const Point* p = boost::get<Point>(&*res)) {
                    intersections.push_back({*p, {i, j}});
                }
            }
        }
    }
    return intersections;
}

double param_along_segment(const Point& a, const Point& b, const Point& p) {
    Vector ab = b - a;
    Vector ap = p - a;
    return (ab * ap) / ab.squared_length(); // projeção escalar normalizada
}


void insert_intersections_and_links(std::vector<Point>& poly, 
                                    const std::vector<std::pair<Point, std::pair<int, int>>>& intersections,
                                    std::vector<JumpLink>& links) {
    std::vector<Point> new_poly;
    std::map<int, std::vector<std::pair<Point, int>>> inserts;

    // Para armazenar diretamente os índices das cópias inseridas
    std::vector<size_t> copy_indices_a(intersections.size());
    std::vector<size_t> copy_indices_b(intersections.size());

    // 1. Marcar quais interseções vão ser inseridas após quais segmentos
    for (size_t k = 0; k < intersections.size(); ++k) {
        const auto& inter = intersections[k];
        inserts[inter.second.first].push_back({inter.first, static_cast<int>(k * 2)});     // cópia 1
        inserts[inter.second.second].push_back({inter.first, static_cast<int>(k * 2 + 1)}); // cópia 2
    }

    // 2. Ordenar as interseções ao longo de cada segmento
    for (auto& [seg_idx, points] : inserts) {
        const Point& a = poly[seg_idx];
        const Point& b = poly[(seg_idx + 1) % poly.size()];
        std::sort(points.begin(), points.end(),
                  [&](const std::pair<Point, int>& lhs, const std::pair<Point, int>& rhs) {
                      return param_along_segment(a, b, lhs.first) < param_along_segment(a, b, rhs.first);
                  });
    }

    // 3. Inserir os pontos no novo polígono, salvando os índices diretamente
    int count = 0;
    for (size_t i = 0; i < poly.size(); ++i) {
        new_poly.push_back(poly[i]);
        count++;

        if (inserts.count(i)) {
            for (auto& ins : inserts[i]) {
                new_poly.push_back(ins.first);

                // Registrar o índice da cópia da interseção k
                int k = ins.second / 2;
                if (ins.second % 2 == 0)
                    copy_indices_a[k] = count;
                else
                    copy_indices_b[k] = count;

                count++;
            }
        }
    }

    // 4. Criar links diretamente dos índices armazenados
    std::set<size_t> used_indices;

    for (size_t k = 0; k < intersections.size(); ++k) {
        size_t idx_a = copy_indices_a[k];
        size_t idx_b = copy_indices_b[k];

        const Point& point_a = new_poly[idx_a];
        const Point& point_b = new_poly[idx_b];

        if (point_a == point_b &&
            used_indices.find(idx_a) == used_indices.end() &&
            used_indices.find(idx_b) == used_indices.end()) {

            links.push_back({idx_a, idx_b});
            used_indices.insert(idx_a);
            used_indices.insert(idx_b);
        } else {
            std::cerr << "⚠️ Problema ao criar link " << k << ": pontos diferentes ou já usados.\n";
        }
    }

    // 5. Mostrar os links criados
    for (size_t k = 0; k < links.size(); ++k) {
        size_t idx_a = links[k].index_a;
        size_t idx_b = links[k].index_b;
        std::cout << "Link " << k << ": idx_a = " << idx_a << " (" 
                  << new_poly[idx_a].x() << ", " << new_poly[idx_a].y() << ") <--> idx_b = " << idx_b << " (" 
                  << new_poly[idx_b].x() << ", " << new_poly[idx_b].y() << ")\n";
    }

    // 6. Substituir o polígono original
    poly = new_poly;
}



std::vector<std::vector<Point>> split_polygon(const std::vector<Point>& poly, const std::vector<JumpLink>& links) {
    std::vector<std::vector<Point>> fragments;
    std::set<int> visited;
    std::map<int, int> jump_map;

    for (const auto& link : links) {
        jump_map[link.index_a] = link.index_b;
        jump_map[link.index_b] = link.index_a;
    }

    int n = poly.size();
    for (int i = 0; i < n; ++i) {
        if (visited.count(i)) continue;

        std::vector<Point> frag;
        int current = i;

        while (!visited.count(current)) {
            frag.push_back(poly[current]);
            visited.insert(current);

            if (jump_map.count(current)) {
                int jump_to = jump_map[current];
                current = (jump_to + 1) % n;
            } else {
                current = (current + 1) % n;
            }
        }

        fragments.push_back(frag);
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

std::vector<Line> offset_polygon(const Polygon& poly, const std::vector<double>& distances) {
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

void export_svg(const Polygon& original, 
                const std::vector<std::vector<Point>>& sub_polygons, 
                const std::string& filename) {
    std::ofstream out(filename);
    out << "<svg xmlns='http://www.w3.org/2000/svg' width='600' height='600' viewBox='-10 -10 120 120'>\n";


    out << "<polygon points='";
    for (auto& p : original)
        out << p.x() << "," << -p.y() << " ";
    out << "' fill='none' stroke='blue' stroke-width='1'/>\n";

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
    Polygon poly;
    poly.push_back(Point(0, 10));
    poly.push_back(Point(5, -10));
    poly.push_back(Point(-5, -5));
    poly.push_back(Point(-5, -20));
    poly.push_back(Point(5, -15));
    poly.push_back(Point(0, -40));
    poly.push_back(Point(20, -40));
    poly.push_back(Point(15, -15));
    poly.push_back(Point(35, -20));
    poly.push_back(Point(35, -5));
    poly.push_back(Point(15, -10));
    poly.push_back(Point(20, 10));
    std::vector<double> distances = {-4.0, -4.0, -4.0, -4.0,-4.0, -4.0, -4.0, -4.0,-4.0, -4.0, -4.0, -4.0};


    auto offset_lines = offset_polygon(poly, distances);
    auto intersections = compute_edge_connections(offset_lines);

    // Debug: imprimir as distâncias dos pontos deslocados
    std::cout << "Distâncias dos deslocamentos:\n";
    for (size_t i = 0; i < offset_lines.size(); ++i) {
        std::cout << "Aresta " << i << ": distância = " << distances[i] << "\n";
    }

    // Debug: imprimir as localizações das interseções
    std::cout << "Interseções do offset:\n";
    for (size_t i = 0; i < intersections.size(); ++i) {
        std::cout << "Interseção " << i << ": (" << intersections[i].x() << ", " << intersections[i].y() << ")\n";
    }

    // Aqui!
    std::cout << "\nPonto original -> Ponto offsetado:\n";
    for (size_t i = 0; i < poly.size(); ++i) {
        if (i < intersections.size()) {
            std::cout << "  (" << poly[i].x() << ", " << poly[i].y() << ")"
                      << "  ->  (" << intersections[i].x() << ", " << intersections[i].y() << ")\n";
        }
    }

    auto internal_inters = find_self_intersections(intersections);

    // Debug: imprimir as interseções internas
    std::cout << "Self-intersections:\n";
    for (size_t i = 0; i < internal_inters.size(); ++i) {
        auto& p = internal_inters[i].first;
        auto& segs = internal_inters[i].second;
        std::cout << "Self-intersection " << i << ": (" << p.x() << ", " << p.y() << ") entre arestas " << segs.first << " e " << segs.second << "\n";
    }

    std::vector<JumpLink> links;
    insert_intersections_and_links(intersections, internal_inters, links);

    // Debug: imprimir a lista dos pontos do polígono com interseções
    std::cout << "Polígono com interseções inseridas:\n";
    for (size_t i = 0; i < intersections.size(); ++i) {
        std::cout << i << ": (" << intersections[i].x() << ", " << intersections[i].y() << ")\n";
    }

    auto sub_polygons = split_polygon(intersections, links);
    auto filtered_polygons = remove_clockwise_polygons(sub_polygons);


    // Debug: imprimir os vértices de cada sub-polígono
    std::cout << "Vértices dos sub-polígonos:\n";
    for (size_t i = 0; i < sub_polygons.size(); ++i) {
        std::cout << "Sub-polígono " << i << ":\n";
        for (size_t j = 0; j < sub_polygons[i].size(); ++j) {
            const auto& p = sub_polygons[i][j];
            std::cout << "  (" << p.x() << ", " << p.y() << ")\n";
        }
    }

    export_svg(poly, filtered_polygons, "output.svg");

    std::cout << "Arquivo 'output.svg' gerado.\n";



    return 0;
}
