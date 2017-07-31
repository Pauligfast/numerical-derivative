//Pavlo Khryshcheniuk
#include <iostream>
#include <cmath>

#include "funkcja.h"

#ifdef sin
#undef sin
#endif
#ifdef cos
#undef cos
#endif
#ifdef exp
#undef exp
#endif

using namespace std;

int k;

class Vector {
public:
    double *x;
public:


    Vector(double x0) {
        x = new double[k];

        x[0] = x0;
        x[1] = 1;
        for (int i = 2; i < k; i++) {
            x[i] = 0;
        }

    }

    Vector() {
        x = new double[k];
        for (int i = 0; i < k; i++) {
            x[i] = 0;
        }

    }

    Vector(const Vector &v) {
        x = new double[k];
        for (int i = 0; i < k; i++) {
            x[i] = v.x[i];
        }

    }

    ~Vector() {

        delete[]x;
    }

    Vector operator+(const Vector &v) const {
        Vector temp(*this);
        for (int i = 0; i < k; i++) {
            temp.x[i] = this->x[i] + v.x[i];
        }
        return temp;
    }


    Vector operator-(const Vector &v) const {
        Vector temp(*this);
        for (int i = 0; i < k; i++) {
            temp.x[i] = this->x[i] - v.x[i];
        }
        return temp;
    }

    Vector operator*(const Vector &v) const {
        Vector temp;
        for (int i = 0; i < k; i++) {
            for (int j = 0; j <= i; j++) {
                temp.x[i] += this->x[i - j] * v.x[j];
            }
        }
        return temp;
    }

    Vector operator=(const Vector &v) {

        for (int i = 0; i < k; i++) {
            this->x[i] = v.x[i];
        }
        return *this;
    }


    Vector operator/(const Vector &v) {
        Vector result;
        double v0 = v.x[0];
        double sum = 0.0;
        for (int i = 0; i < k; i++) {
            sum = 0.0;

            for (int j = 1; j <= i; j++) {
                sum += v.x[j] * result.x[i - j];
            }
            result.x[i] = ((1.0 / v0) * (this->x[i] - sum));
        }
        return result;
    }

    Vector operator-() const {
        Vector temp;
        for (int i = 0; i < k; i++) {
            temp.x[i] = -(this->x[i]);
        }
        return temp;
    }


    friend ostream &operator<<(ostream &os, const Vector &v);

    friend Vector operator+(const double left, const Vector &right);

    friend Vector operator+(const Vector &left, const double right);

    friend Vector operator-(const double left, const Vector &right);

    friend Vector operator-(const Vector &left, const double right);

    friend Vector operator*(const double left, const Vector &right);

    friend Vector operator*(const Vector &left, const double right);

    friend Vector operator/(const double left, const Vector &right);

    friend Vector operator/(const Vector &left, const double right);


};

Vector sin(const Vector &v) {
    Vector temp1, temp2;
    temp1.x[0] = sin(v.x[0]);
    temp2.x[0] = cos(v.x[0]);
    for (int i = 1; i < k; i++) {
        for (int j = 0; j < i; j++) {
            temp1.x[i] += (i - j) * temp2.x[j] * v.x[i - j];
        }
        temp1.x[i] /= i;
        for (int j = 0; j < i; j++) {
            temp2.x[i] += (i - j) * temp1.x[j] * v.x[i - j];
        }
        temp2.x[i] /= -i;
    }
    return temp1;
}

Vector cos(const Vector &v) {
    Vector temp1, temp2;
    temp1.x[0] = cos(v.x[0]);
    temp2.x[0] = sin(v.x[0]);
    for (int i = 1; i < k; i++) {
        for (int j = 0; j < i; j++) {
            temp1.x[i] += (i - j) * temp2.x[j] * v.x[i - j];
        }
        temp1.x[i] /= -i;
        for (int j = 0; j < i; j++) {
            temp2.x[i] += (i - j) * temp1.x[j] * v.x[i - j];
        }
        temp2.x[i] /= i;
    }
    return temp1;
}

Vector exp(const Vector &v) {
    Vector temp;
    temp.x[0] = exp(v.x[0]);
    for (int i = 1; i < k; i++) {
        for (int j = 0; j < i; j++) {
            temp.x[i] += (i - j) * temp.x[j] * v.x[i - j];

        }
        temp.x[i] /= i;
    }
    return temp;
}

ostream &operator<<(ostream &os, const Vector &v) {
    long factorial = 1;
    cout.precision(16);
    cout << fixed << v.x[0] << " ";
    for (int i = 1; i < k - 1; i++) {
        factorial *= i;
        cout << v.x[i] * factorial << " ";

    }
    factorial *= k - 1;
    cout << v.x[k - 1] * factorial << endl;
}

Vector operator+(const double left, const Vector &right) {
    Vector temp(left);
    temp.x[1] = 0;
    for (int i = 0; i < k; i++) {
        temp.x[i] += right.x[i];
    }
    return temp;
}

Vector operator+(const Vector &left, const double right) {
    Vector temp(right);
    temp.x[1] = 0;

    for (int i = 0; i < k; i++) {
        temp.x[i] += left.x[i];
    }
    return temp;
}

Vector operator-(const double left, const Vector &right) {
    Vector temp(left);
    temp.x[1] = 0;

    for (int i = 0; i < k; i++) {
        temp.x[i] -= right.x[i];
    }
    return temp;
}

Vector operator-(const Vector &left, const double right) {
    Vector temp(right);
    Vector lefttemp(left);
    temp.x[1] = 0;
    for (int i = 0; i < k; i++) {
        lefttemp.x[i] = left.x[i] - temp.x[i];
    }
    return lefttemp;
}

Vector operator*(const Vector &left, const double right) {
    Vector temp;
    Vector rightV(right);
    rightV.x[1] = 0;
    for (int i = 0; i < k; i++) {
        for (int j = 0; j <= i; j++) {
            temp.x[i] += left.x[i - j] * rightV.x[j];
        }
    }
    return temp;
}

Vector operator*(const double left, const Vector &right) {
    Vector temp;
    Vector leftV(left);
    leftV.x[1] = 0;
    for (int i = 0; i < k; i++) {
        for (int j = 0; j <= i; j++) {
            temp.x[i] += leftV.x[i - j] * right.x[j];
        }
    }
    return temp;
}

/*Vector operator/(const Vector &left, const double right) {
    Vector rightV(right);
    rightV.x[1] = 0;
    Vector temp;
    for (int i = 0; i < k; i++) {
        for (int j = 1; j <= i; j++) {
            temp.x[i] += rightV.x[j] * (left.x[i - j] / rightV.x[i - j]);
        }
        temp.x[i] = (left.x[i] - temp.x[i]) / rightV.x[0];
    }
    return temp;
}*/

Vector operator/(const Vector &v, const double right) {
    Vector r(right);
    r.x[1] = 0;
    Vector temp;
    temp.x[0] = v.x[0] / r.x[0];
    for (int i = 1; i < k; ++i) {
        temp.x[i] = v.x[i];
        for (int j = 1; j <= i; ++j) {
            temp.x[i] -= r.x[j] * temp.x[i - j];
        }
        temp.x[i] /= r.x[0];
    }
    return temp;
}

Vector operator/(const double left, const Vector &r) {
    Vector v(left);
    r.x[1] = 0;
    Vector temp;
    temp.x[0] = v.x[0] / r.x[0];
    for (int i = 1; i < k; ++i) {
        temp.x[i] = v.x[i];
        for (int j = 1; j <= i; ++j) {
            temp.x[i] -= r.x[j] * temp.x[i - j];
        }
        temp.x[i] /= r.x[0];
    }
    return temp;
}

/*Vector operator/(const double left, const Vector &right) {
    Vector leftV(left);
    leftV.x[1] = 0;
    Vector temp;
    for (int i = 0; i < k; i++) {
        for (int j = 1; j <= i; j++) {
            temp.x[i] += right.x[j] * (leftV.x[i - j] / right.x[i - j]);
        }
        temp.x[i] = (leftV.x[i] - temp.x[i]) / right.x[0];
    }
    return temp;
}*/


int main() {


    int N, M;

    cin >> N;

    cin >> M;
    double x[M];
    k = N + 1;
    for (int i = 0; i < M; i++) {

        cin >> x[i];

    }
    for (int i = 0; i < M; i++) {
        Vector v(x[i]);
        v = funkcja(v);
        cout << v;

    }
    // cin >> N;
    return 0;
}
