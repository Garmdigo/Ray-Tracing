// C++ include
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
//https://stackoverflow.com/questions/1727881/how-to-use-the-pi-constant-in-c
//To get PI
//https://stackoverflow.com/questions/8690567/setting-an-int-to-infinity-in-c
// Making an infinite variable
#include <limits>
#include <math.h>
#define _USE_MATH_DEFINES
// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"
#include "utils.h"

// Shortcut to avoid Eigen:: and std:: everywhere, DO NOT USE IN .h
using namespace std;
using namespace Eigen;
//http://paulbourke.net/geometry/polygonmesh/
//http://what-when-how.com/advanced-methods-in-computer-graphics/mesh-processing-advanced-methods-in-computer-graphics-part-1/
//https://www.youtube.com/watch?v=P2xMqTgsgsE


//class Mesh
//{
//public:
//	Mesh(string f)
//	{
//		string fileName;
//		int temp;
//		Vector3f zero(0, 0, 0);
//		ifstream F(f);
//		F >> fileName >> numVertices >> numTriangles >> temp;
//
//		for (int i = 0; i < numVertices; i++)
//		{
//			float x, y, z;
//			F >> x >> y >> z;
//			Vector3f temp(x, y, z);
//			V.push_back(temp);
//			NAngle.push_back(zero);
//			NA.push_back(zero);
//		}
//
//
//		for (int i = 0; i < numTriangles; i++)
//		{
//			int temp2, temp3, temp4;
//			F >> temp >> temp2 >> temp3 >> temp4;
//			Vector3i temp(temp2, temp3, temp4);
//			T.push_back(temp);
//		}
//
//		for (int i = 0; i < numTriangles; i++)
//		{
//			auto x = V[T[i].x()];
//			auto y = V[T[i].y()];
//			auto z = V[T[i].z()];
//			Vector3f temp = x;
//			Vector3f temp2 = y;
//			Vector3f temp3 = z;
//
//			Vector3f temp4 = temp3 - temp;
//			Vector3f temp5 = temp2 - temp;
//
//			Vector3f N = temp5.cross(temp4);
//			normal.push_back(N);
//		}
//
//		for (int i = 0; i < numTriangles; i++)
//		{
//			NA[T[i].x()] += normal[i];
//			NA[T[i].y()] += normal[i];
//			NA[T[i].z()] += normal[i];
//		}
//
//		for (int i = 0; i < numTriangles; i++)
//		{
//			normal[i] = (normal[i] / normal[i].norm());
//		}
//		for (int i = 0; i < numVertices; i++)
//		{
//			NA[i] = NA[i] / NA[i].norm();
//		}
//

//
//

//				
//
//			}
//
//	}
//	vector<Vector3f> getV()
//	{
//		return V;
//	}
//	vector<Vector3f> getN()
//	{
//		return normal;
//	}
//	vector<Vector3f> getNArea()
//	{
//		return NA;
//	}
//	vector<Vector3f> getNAngle()
//	{
//		return NAngle;
//	}
//
//	vector<Vector3i> getT()
//	{
//		return T;
//	}
//	int getNumT()
//	{
//		return numTriangles;
//	}
//	int getNumV()
//	{
//		return numVertices;
//	}
//	vector<Vector3f> getBoundingBox()
//	{
//		return bounding_box;
//	}
//
//private:
//	int numVertices;
//	int numTriangles;
//	vector<Vector3f> V;
//	vector<Vector3f > normal;
//	vector<Vector3f > NA; //Normal Area
//	vector<Vector3f>  NAngle;
//	vector<Vector3i>  T;
//	vector<Vector3f>  bounding_box;
//
//};




class Illumination
{
public:
	Illumination(const Vector3d &p, const double(&i))
	{
		P = p;
		I = i;
	}
	Vector3d getP() const
	{
		return P;
	}
	double getI() const
	{
		return I;
	}
private:
	Vector3d P;
	double I;
};

class Ray
{
public:
	Ray(const Vector3d & O, const Vector3d  & dir)
	{
		Origin = O, Direction = dir;
	}
	auto getOrigin() const
	{
		return Origin;
	}
	auto getDirection() const
	{
		return Direction;
	}
private:
	Vector3d Origin;
	Vector3d Direction;
};

class Description
{
public:
	Description(const double &r, const Vector3d &al, const double &s)
	{

		Albedo = al;
		Refract = r;
		Spec = s;
	}
	Description()
	{
		Refract = 1;
		Albedo = Vector3d(1, 0, 0);

	}

	auto getAlbedo() const
	{
		return Albedo;
	}
	auto getRefract() const
	{
		return Refract;
	}

	auto getSpec() const
	{
		return Spec;
	}
	auto setAlbedo(Vector3d b)
	{
		Albedo = b;
	}
private:
	double Refract;
	Vector3d Albedo;
	double Spec;


};

class Sphere
{
public:

	Sphere(const Vector3d &c, const float &r, const Description &D)
	{
		Object = D;
		center = c;
		radius = r;
	}

	bool intersect(Ray R)
	{
		auto O = R.getOrigin() - center;
		auto D = R.getDirection();
		double DO = D.dot(O);
		double OO = O.dot(O);
		double DD = D.dot(D);

		if (DD <= 0) {
			return false;
		}

		double Delta = 4 * DO * DO - 4 * DD * (OO - radius * radius);
		double t1 = (-2 * DO - sqrt(Delta)) / (2 * DD);
		double t2 = (-2 * DO + sqrt(Delta)) / (2 * DD);

		if (t1 < t2) {
			t = t1;
		}
		else {
			t = t2;
		}
		return true;
	}
	auto getCenter()
	{
		return center;
	}
	auto getRadius()
	{
		return radius;
	}
	auto getT()
	{
		return t;
	}
	auto getDescription()
	{
		return Object;
	}

	double t = 0;
private:
	Vector3d center;
	double radius;
	Description Object;
};
class Plane {
public:
	Plane() : Origin{ 0, -1, 0 }, Normal{ 1, 0, 0 } {}


	bool getI(const Ray &ray)
	{
		double denom = Normal.dot(ray.getDirection());
		if (fabs(denom) > 0) {
			auto temp = (Origin - ray.getOrigin()).dot(Normal);
			double temp2 = temp / denom;
			if (temp2 > 0)
			{
				return (t >= 0);
			}
		}
		return false;
	}
	Vector3d GetOrigin() const
	{
		return Origin;
	}

	Vector3d getN(const Vector3d &point)
	{
		return Normal;
	}
private:
	Vector3d Normal, Origin;
	double t;
};
//I−2(N⋅I)N
Vector3d Reflection(const Vector3d &I, const Vector3d &N)
{
	return I - 2 * N.dot(I) * N;
}
bool SphereIntersect(Ray & ray, Vector3d &hit, Vector3d &N, vector<Sphere> &S, Description &Des)
{
	Ray R = ray;


	auto minT = numeric_limits<double>::infinity();
	int minI = -1;
	for (int i = 0; i < S.size(); i++)
	{
		if (S[i].intersect(R) && S[i].getT() < minT)
		{
			minT = S[i].getT();
			minI = i;
		}
	}
	if (minI < 0) {
		return false;
	}
	else {
		hit = R.getOrigin() + R.getDirection() * minT;
		N = (hit - S[minI].getCenter()).normalized();
		Des = S[minI].getDescription();

	}

	return true;
}



Vector3d Shade(const Description& material, const Illumination& light, Vector3d &normal, Vector3d & hit, Vector3d & origin)
{
	const auto e = 1e-3;
	auto l = (light.getP() - hit).normalized();
	auto lDistance = (light.getP() - hit).normalized();


	auto cosAngle = max(l.dot(normal), 0.0);
	auto diffuse = material.getAlbedo() * cosAngle;

	auto v = (origin - hit).normalized();
	auto rv = Reflection(-v, normal);
	cosAngle = max(l.dot(v), 0.0);
	auto specular = pow(cosAngle, material.getSpec()) * Vector3d(1, 1, 1);

	return (diffuse + specular) * light.getI();
}

Vector3d Shade(const Description& material, const vector<Illumination>& lights, Vector3d normal, Vector3d hit, Vector3d origin) {
	Vector3d sum = Vector3d(0, 0, 0);
	Vector3d point, N;
	for (auto l : lights) {

		sum += Shade(material, l, normal, hit, origin);
	}
	return sum;
}

void Tracing(const vector <Sphere> &S, vector<Illumination> & Light, Plane Plane, bool useOrthoProjection, string Filename)
{
	//	Mesh Bunny("bunny.off");
	//	Mesh Cube("bumpy_cube.off");

	cout << "Running Ray tracing program" << endl;
	auto P = Plane;
	auto Sphere = S;

	double Width = 640;
	double Height = 640;
	double Dimension = Width * Height;
	vector <Vector3d> Display(Dimension);

	/*bool useOrthoProjection = true;*/

	for (unsigned W = 0; W < Width; W++)
	{
		for (unsigned L = 0; L < Height; L++)
		{
			// Convert (W,L) to (-1, 1)
			double X = (L + 0.5) / (Height / 2) - 1;
			double Y = (W + 0.5) / (Width / 2) - 1;
			const double focalLength = 0.9;

			const Vector3d origin = useOrthoProjection ? Vector3d(X, Y, focalLength) : Vector3d(0, 0, focalLength);
			const Vector3d pixelPosition = Vector3d(X, Y, 0);

			auto rayDirection = (pixelPosition - origin).normalized();

			Ray R(origin, rayDirection);
			Vector3d point = Vector3d(.07, .2, .1);
			Vector3d N;
			Description D;
			
			if (P.getI(R)>0)
			{
				Display[L + W * Width] = Vector3d(0, 0, 1);
			}

			if (SphereIntersect(R, point, N, Sphere, D))
			{
				Display[L + W * Width] = Shade(D, Light, N, point, origin);
			}
			else
			{
				Display[L + W * Width] = Vector3d(40, 40, 40);
			}
		}
	}

	MatrixXd R = MatrixXd::Zero(Height, Width); // Store A
	MatrixXd G = MatrixXd::Zero(Height, Width); // Store G
	MatrixXd B = MatrixXd::Zero(Height, Width); // Store B
	MatrixXd A = MatrixXd::Zero(Height, Width); // Store A

	for (unsigned W = 0; W < Width; W++)
	{
		for (unsigned L = 0; L < Height; L++)
		{
			R(W, L) = Display[W + L * Width].x();
			G(W, L) = Display[W + L * Width].y();
			B(W, L) = Display[W + L * Width].z();
			A(W, L) = 1;
		}
	}

	write_matrix_to_png(R, G, B, A, Filename);
}

int main()
{
	Plane P;
	vector<Sphere> S;
	Description Object1(1.0, Vector3d(0.6, 0.3, 0.1), 100);
	S.push_back(Sphere(Vector3d(0, 0, -6), 1, Object1));

	Description Object2(1.0, Vector3d(0.1, 0.9, 0.1), 100);
	S.push_back(Sphere(Vector3d(-3, 2, -4), 1, Object2));

	Description Object3(1.0, Vector3d(0.3, 0.2, 0.7), 1000);
	S.push_back(Sphere(Vector3d(4, -3, -5), 1, Object3));

	Description Object4(1.0, Vector3d(0.3, 0.2, 0.7), 1000);
	S.push_back(Sphere(Vector3d(0, -200000, -5), 1, Object4));

	vector<Illumination> L;
	L.push_back(Illumination(Vector3d(1, 0, 5), 0.8));
	L.push_back(Illumination(Vector3d(5, 5, 5), 0.7));
	L.push_back(Illumination(Vector3d(15, 15, -5), 0.6));

	Tracing(S, L,P, 1, "OrthoProjection.png");
	Tracing(S, L,P,0, "PerspectiveProjection.png");
	cout << "Program is over! Goodbye!";
	return 0;
}

