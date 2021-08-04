
typedef struct {
    double a;
    double b;
    double c;
} vector_t;

void vector_init(vector_t *v, double a, double b, double c);
void vector_add(vector_t *v, vector_t *add);
double vector_get_magnitude(vector_t *v);
void vector_set_magnitude(vector_t *v, double new_magnitude);


