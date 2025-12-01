#include <iostream>
#include <vector>
#include <random>
#include <chrono>

//Gerador de grafo aleatório
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>

using namespace boost;
using namespace std;

//Grafo com peso nas arestas 
using Grafo = adjacency_list<vecS, vecS, undirectedS,
                             no_property,
                             property<edge_weight_t, double>>;

using Vertice = graph_traits<Grafo>::vertex_descriptor;
// Gera um grafo aleatório com n vértices e densidade (0 a 1). ja atribuindo peso ao criar a aresta)
Grafo gerarGrafo(int n, double densidade, double wmin = 1.0, double wmax = 10.0) {
    Grafo g(n);

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> prob(0.0, 1.0);
    uniform_real_distribution<> weight(wmin, wmax);

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (prob(gen) < densidade) {
                double w = weight(gen);
                add_edge(i, j, w, g); // já definimos o peso
            }
        }
    }

    return g;
}


// Roda Dijkstra e retorna tempo em milissegundos.
// Também preenche vetores de distância e predecessor
double rodarDjikstra(const Grafo& g, int start, vector<double>* out_dist = nullptr, vector<int>* out_pred = nullptr) {
    int n = static_cast<int>(num_vertices(g));
    vector<double> dist(n, (numeric_limits<double>::max)());
    vector<int> pred(n, -1);

    // criar property maps a partir de vetores
    auto index_map = get(vertex_index, g);
    auto dist_map = make_iterator_property_map(dist.begin(), index_map);
    auto pred_map = make_iterator_property_map(pred.begin(), index_map);

    auto t0 = chrono::high_resolution_clock::now();

    dijkstra_shortest_paths(
        g,
        start,
        predecessor_map(pred_map).distance_map(dist_map)
    );

    auto t1 = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> ms = t1 - t0;

    if (out_dist) *out_dist = std::move(dist);
    if (out_pred) *out_pred = std::move(pred);

    return ms.count();
}

// função principal-

int main() {
    int n;
    double densidade;

    cout << "Numero de vertices: ";
    if (!(cin >> n)) return 0;
    cout << "Densidade (0 a 1): ";
    if (!(cin >> densidade)) return 0;

    cout << "Gerando grafo (" << n << " vertices, densidade=" << densidade << ")...\n";
    Grafo g = gerarGrafo(n, densidade);

    cout << "Rodando Dijkstra a partir do vertice 0...\n";
    vector<double> dist;
    vector<int> pred;
    double tempo = rodarDjikstra(g, 0, &dist, &pred);

    cout << "Tempo: " << tempo << " ms\n\n" << endl;

    cout << "Distancias (ate " << n << " vertices):\n";
    for (int i = 0; i < n; ++i) {
        if (dist[i] == (numeric_limits<double>::max)())
            cout << i << ": INF\n";
        else
            cout << i << ": " << dist[i] << "\n";
    }

    return 0;
}