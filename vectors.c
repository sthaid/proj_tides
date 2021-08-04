#include <vectors.h>
#include <math.h>
#include <assert.h>

static inline double square(double x) { return x * x; }

void vector_init(vector_t *v, double a, double b, double c)
{
    v->a = a;
    v->b = b;
    v->c = c;
}

void vector_add(vector_t *v, vector_t *add)
{
    v->a += add->a;
    v->b += add->b;
    v->c += add->c;
}

double vector_get_magnitude(vector_t *v)
{
    return sqrt( square(v->a) + square(v->b) + square(v->c) );
}

void vector_set_magnitude(vector_t *v, double new_magnitude)
{
    double current_magnitude, factor;

    current_magnitude = vector_get_magnitude(v);
    assert(current_magnitude > 0);

    factor = new_magnitude / current_magnitude;

    v->a *= factor;
    v->b *= factor;
    v->c *= factor;
}
