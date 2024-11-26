#include "Matrix.h"
#include "Transform.h"
class Pivot {
    public:
    Pivot(vector<Transform> FrameTransformationList);
    Matrix p_t;
    Matrix p_post;
};