// I have memory leaks where I use `new` 
// I use very inefficient algorithm for finding neighbours of tiles (naive approach)

#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <math.h>
#include <map>

#define PI 3.14159265358979323846

using namespace std;

typedef pair<float, float> Point;

const int N_ITERATIONS = 4; // Number of recursion iterations

const vector<float> IDENTITY {1, 0, 0, 0, 1, 0};
const vector<string> TILE_NAMES {"Gamma", "Delta", "Theta", "Lambda", "Xi", "Pi", "Sigma", "Phi", "Psi"};
const vector<Point> SPECTRE_POINTS {
    Point {0, 0},
    Point {1, 0},
    Point {1.5, -sqrt(3)/2},
    Point {1.5+sqrt(3)/2, 0.5-sqrt(3)/2},
    Point {1.5+sqrt(3)/2, 1.5-sqrt(3)/2},
    Point {2.5+sqrt(3)/2, 1.5-sqrt(3)/2},
    Point {3+sqrt(3)/2,   1.5},
    Point {3, 2},
    Point {3-sqrt(3)/2,   1.5},
    Point {2.5-sqrt(3)/2, 1.5+sqrt(3)/2},
    Point {1.5-sqrt(3)/2, 1.5+sqrt(3)/2},
    Point {0.5-sqrt(3)/2, 1.5+sqrt(3)/2},
    Point {-sqrt(3)/2,    1.5},
    Point {0, 1}
};

int num_tiles = 0;
vector<vector<Point>> tiles_vertices;

// Affine matrix multiply
vector<float> mul(const vector<float> &A, const vector<float> &B) {
    vector<float> C {
        A[0]*B[0] + A[1]*B[3],
        A[0]*B[1] + A[1]*B[4],
        A[0]*B[2] + A[1]*B[5] + A[2],

        A[3]*B[0] + A[4]*B[3],
        A[3]*B[1] + A[4]*B[4],
        A[3]*B[2] + A[4]*B[5] + A[5]
    };
    return C;
}

// Rotation matrix
vector<float> trot(float ang) {
    float c = cos(ang);
    float s = sin(ang);
    vector<float> M {c, -s, 0, s, c, 0};
    return M;
}

// Translation matrix
vector<float> ttrans(float tx, float ty) {
    vector<float> M {1, 0, tx, 0, 1, ty};
    return M;
}

vector<float> transTo(const Point &p, const Point &q) {
    return ttrans(q.first - p.first, q.second - p.second);
}

// Matrix * point
Point transPt(const vector<float> &M, const Point &P) {
    Point newP {M[0]*P.first + M[1]*P.second + M[2], M[3]*P.first + 
        M[4]*P.second + M[5]};
    return newP;
}


class BaseTile {
protected:
    vector<Point> quad;
    vector<float> transform;
public:
    virtual void setTransform(const vector<float> &transf) = 0;
    const vector<Point> getQuad() {
        return this->quad;
    }
};


class Tile: public BaseTile {
private:
    vector<Point> points;

public:
    Tile(const vector<Point> &points) {
        //points: list of Tile coordinate points
        this->quad = {points[3], points[5], points[7], points[11]};
        this->points = points;
    }

    void setTransform(const vector<float> &transf) {
        this->transform = transf;

        vector<Point> finalPoints(this->points.size());
        for (int i = 0; i < this->points.size(); ++i) {
            finalPoints[i] = transPt(transf, this->points[i]);
        }
        tiles_vertices.push_back(finalPoints);
    }
};


class MetaTile: public BaseTile {
private:
    vector<pair<BaseTile*, vector<float>>> geometries;

public:
    MetaTile(const vector<pair<BaseTile*, vector<float>>> &geometries,
                const vector<Point> &quad) {
        /*
        geometries: list of pairs of (Meta)Tiles and their transformations
        quad: MetaTile quad points
        */
        this->geometries = geometries;
        this->quad = quad;
    }

    void setTransform(const vector<float> &transf) {
        this->transform = transf;
        for (pair<BaseTile*, vector<float>> tileAndTransf: this->geometries) {
            BaseTile* tile = tileAndTransf.first;
            vector<float> tileTransf = tileAndTransf.second;
            tile->setTransform(mul(transf, tileTransf));
        }
    }
};


map<string, BaseTile*> buildSpectreBase() {
    map<string, BaseTile*> spectre_base_cluster;
    for (string tileName : TILE_NAMES) {
        if (tileName != "Gamma") {
            spectre_base_cluster[tileName] = new Tile(SPECTRE_POINTS);
        }
    }
    // special rule for Gamma
    BaseTile* gamma1 = new Tile(SPECTRE_POINTS);
    BaseTile* gamma2 = new Tile(SPECTRE_POINTS);
    pair<BaseTile*, vector<float>> gamma1Geom {gamma1, IDENTITY};
    vector<float> gamma2Transf = mul(ttrans(SPECTRE_POINTS[8].first,
        SPECTRE_POINTS[8].second), trot(PI/6));
    pair<BaseTile*, vector<float>> gamma2Geom {gamma2, gamma2Transf};
    vector<pair<BaseTile*, vector<float>>> geoms {gamma1Geom, gamma2Geom};
    vector<Point> gamma_quad {SPECTRE_POINTS[3], SPECTRE_POINTS[5], SPECTRE_POINTS[7], SPECTRE_POINTS[11]};
    MetaTile* mystic = new MetaTile(geoms, gamma_quad);
    
    spectre_base_cluster["Gamma"] = mystic;
    return spectre_base_cluster;
}

map<string, BaseTile*> buildSupertiles(map<string, BaseTile*> &tileSystem) {
    /*   
    iteratively build on current system of tiles
    tileSystem = current system of tiles, initially built with buildSpectreBase()
    */

    // First, use any of the nine-unit tiles in tileSystem to obtain
    // a list of transformation matrices for placing tiles within
    // supertiles.
    vector<Point> quad = tileSystem["Delta"]->getQuad();
    vector<float> R {-1, 0, 0, 0, 1, 0};

    //[rotation angle, starting quad point, target quad point]
    vector<vector<float>> transformation_rules = {
        {60, 3, 1}, {0, 2, 0}, {60, 3, 1}, {60, 3, 1},
        {0, 2, 0}, {60, 3, 1}, {-120, 3, 3}
    };

    vector<vector<float>> transformations {IDENTITY};
    float total_angle = 0;
    vector<float> rotation = IDENTITY;
    vector<Point> transformed_quad = quad;

    for (vector<float> transf_rule: transformation_rules) {
        float _angle = transf_rule[0];
        float _from = transf_rule[1];
        float _to = transf_rule[2];

        if (_angle != 0) {
            total_angle += _angle;
            rotation = trot(total_angle * PI / 180);
            for (int i = 0; i < 4; ++i) {
                transformed_quad[i] = transPt(rotation, quad[i]);
            }
        }

        vector<float> ttt = transTo(
            transformed_quad[_to],
            transPt(transformations.back(), quad[_from])
        );
        transformations.push_back(mul(ttt, rotation));
    }

    for (int i = 0; i < transformations.size(); ++i) {
        transformations[i] = mul(R, transformations[i]);
    }
    
    // Now build the actual supertiles, labelling appropriately.
    map<string, vector<string>> super_rules {
        {"Gamma",  {"Pi",  "Delta", "",  "Theta", "Sigma", "Xi",  "Phi",    "Gamma"}},
        {"Delta",  {"Xi",  "Delta", "Xi",  "Phi",   "Sigma", "Pi",  "Phi",    "Gamma"}},
        {"Theta",  {"Psi", "Delta", "Pi",  "Phi",   "Sigma", "Pi",  "Phi",    "Gamma"}},
        {"Lambda", {"Psi", "Delta", "Xi",  "Phi",   "Sigma", "Pi",  "Phi",    "Gamma"}},
        {"Xi",     {"Psi", "Delta", "Pi",  "Phi",   "Sigma", "Psi", "Phi",    "Gamma"}},
        {"Pi",     {"Psi", "Delta", "Xi",  "Phi",   "Sigma", "Psi", "Phi",    "Gamma"}},
        {"Sigma",  {"Xi",  "Delta", "Xi",  "Phi",   "Sigma", "Pi",  "Lambda", "Gamma"}},
        {"Phi",    {"Psi", "Delta", "Psi", "Phi",   "Sigma", "Pi",  "Phi",    "Gamma"}},
        {"Psi",    {"Psi", "Delta", "Psi", "Phi",   "Sigma", "Psi", "Phi",    "Gamma"}}
    };
    vector<Point> super_quad = {
        transPt(transformations[6], quad[2]),
        transPt(transformations[5], quad[1]),
        transPt(transformations[3], quad[2]),
        transPt(transformations[0], quad[1])
    };

    map<string, BaseTile*> superTileSystem;
    for (auto const& [label, substitutions]: super_rules) {
        vector<pair<BaseTile*, vector<float>>> geoms;

        for (int i = 0; i < transformations.size(); ++i) {
            string substitution = substitutions[i];
            vector<float> transformation = transformations[i];
            if (substitution != "") {
                pair<BaseTile*, vector<float>> geom {
                    tileSystem[substitution],
                    transformation
                };
                geoms.push_back(geom);
            }
        }

        superTileSystem[label] = new MetaTile(geoms, super_quad);
    }
    return superTileSystem;
}


class TileCell {
private:
    vector<Point> verteces;
    vector<int> neighbours;
    vector<int> vertNeighbours;
public:
    TileCell(const vector<Point> &verteces) {
        this->verteces = verteces;
    }

    void addNeighbour(int neighbour) {
        this->neighbours.push_back(neighbour);
    }

    void addVertNeighbour(int neighbour) {
        this->vertNeighbours.push_back(neighbour);
    }

    void print() {
        for (int i = 0; i < this->verteces.size(); ++i) {
            std::cout << "(" + to_string(this->verteces[i].first) + ", " + 
                        to_string(this->verteces[i].second) + ")";
        }
        std::cout << "\n";
        for (int i = 0; i < this->neighbours.size(); ++i) {
            std::cout << this->neighbours[i] << " ";
        }
        std::cout << "\n";
        for (int i = 0; i < this->vertNeighbours.size(); ++i) {
            std::cout << this->vertNeighbours[i] << " ";
        }
        std::cout << "\n";
    }
};


bool equalPoints(const Point &point1, const Point &point2) {
    float eps = 0.001;
    return ((abs(point1.first - point2.first) < eps) 
        && (abs(point1.second - point2.second) < eps));
}

bool haveCommonEdges(const vector<pair<Point, Point>> &edges1, 
                     const vector<pair<Point, Point>> &edges2) {
    for (pair<Point, Point> edge1: edges1) {
        for (pair<Point, Point> edge2: edges2) {
            if ((equalPoints (edge1.first, edge2.first) && 
                 equalPoints (edge1.second, edge2.second)) ||
                (equalPoints (edge1.first, edge2.second) && 
                 equalPoints (edge1.second, edge2.first))) {
                return true;
            }
        }
    }
    return false;
}

bool haveCommonVerts(const vector<Point> &verts1, 
                     const vector<Point> &verts2) {
    for (Point vert1: verts1) {
        for (Point vert2: verts2) {
            if (equalPoints(vert1, vert2)) {
                return true;
            }
        }
    }
    return false;
} 


int main() {
    map<string, BaseTile*> shapes = buildSpectreBase();
    for (int i = 0; i < N_ITERATIONS; ++i) {
        shapes = buildSupertiles(shapes);
    }
    shapes["Delta"]->setTransform(IDENTITY);

    vector<vector<pair<Point, Point>>> tiles_edges(tiles_vertices.size());
    for (int i = 0; i < tiles_vertices.size(); ++i) {
        vector<Point> tile_vertices = tiles_vertices[i];
        vector<pair<Point, Point>> tile_edges;
        for (int j = 0; j < tile_vertices.size()-1; ++j) {
            tile_edges.push_back(make_pair(tile_vertices[j], tile_vertices[j+1]));
        }
        tile_edges.push_back(make_pair(tile_vertices[tile_vertices.size()-1], tile_vertices[0]));
        tiles_edges[i] = tile_edges;
    }

    vector<TileCell*> tileCells(tiles_vertices.size());
    for (int i = 0; i < tiles_vertices.size(); ++i) {
        tileCells[i] = new TileCell(tiles_vertices[i]);
    }

    for (int i = 0; i < tiles_edges.size(); ++i) {
        vector<pair<Point, Point>> tile1Edges = tiles_edges[i];
        vector<Point> tile1Verts = tiles_vertices[i];
        for (int j = i + 1; j < tiles_edges.size(); ++j) {
            vector<pair<Point, Point>> tile2Edges = tiles_edges[j];
            vector<Point> tile2Verts = tiles_vertices[j];
            if (haveCommonEdges(tile1Edges, tile2Edges)) {
                tileCells[i]->addNeighbour(j);
                tileCells[j]->addNeighbour(i);
                continue;
            }
            if (haveCommonVerts(tile1Verts, tile2Verts)) {
                tileCells[i]->addVertNeighbour(j);
                tileCells[j]->addVertNeighbour(i);
            }
        }
    }

    for (int i = 0; i < tileCells.size(); ++i) {
        tileCells[i]->print();
    }

    return 0;
}
