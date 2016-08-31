#ifndef SPS_H
#define SPS_H
struct CSRGraph {
  int nnodes;
  int nedges;
  int *row_offsets;
  int *columns;
  bool *sat;
  float *bias;
  bool *value;
  
  bool alloc() {
    assert(nnodes > 0);
    assert(nedges > 0);
    row_offsets = (int *) calloc(nnodes + 1, sizeof(*row_offsets));
    columns = (int *) calloc(nedges, sizeof(*columns));
    sat = (bool *) calloc(nnodes, sizeof(bool));
    bias = (float *) calloc(nnodes, sizeof(float));
    value = (bool *) calloc(nnodes, sizeof(bool));
    return (row_offsets != NULL) && (columns != NULL) && sat && bias && value;
  }
  
  void set_last_offset() {
    row_offsets[nnodes] = nedges;
  }

  int degree(const int node) const {
    return row_offsets[node + 1] - row_offsets[node];
  }

  void dump_edges() const {
    int i;
    for(i = 0; i < nedges; i++)
      printf("%d ", columns[i]);
    printf("\n");
  }

  void dump_offsets() const {
    int i;
    for(i = 0; i < nnodes; i++)
      printf("%d ", row_offsets[i]);
    printf("\n");
  }
};

struct Edge {
  int nedges;
  int *src;
  int *dst;
  bool *bar;
  float *eta;
  float *pi_0;
  float *pi_S;
  float *pi_U;

  bool alloc() {
    assert(nedges > 0);
    src = (int*) calloc(nedges, sizeof(*src));
    dst = (int*) calloc(nedges, sizeof(*dst));
    bar = (bool*) calloc(nedges, sizeof(*bar));
    eta = (float*) calloc(nedges, sizeof(*eta));
    pi_0 = (float*) calloc(nedges, sizeof(*pi_0));
    pi_S = (float*) calloc(nedges, sizeof(*pi_S));
    pi_U = (float*) calloc(nedges, sizeof(*pi_U));
    return (src && dst && bar && eta && pi_0 && pi_S && pi_U);
  }
};
#endif
