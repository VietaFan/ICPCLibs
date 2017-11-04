// mostly similar to code from Computational Geometry in C (O'Rourke 1994)

class Point {
	public double x, y;
	public Point() {
		this.x = 0;
		this.y = 0;
	}
	public Point(double x, double y) {
		this.x = x;
		this.y = y;
	}
	public double dist(Point other) {
		return Math.sqrt((this.x-other.x)*(this.x-other.x)+(this.y-other.y)*(this.y-other.y));
	}
}

public class LibGeom {
	// Shoelace Formula, 3 var case
	double triArea(Point A, Point B, Point C) {
		return 0.5*Math.abs(A.x*B.y+B.x*C.y+C.x*A.y-A.y*B.x+B.y*C.x+C.y*A.x);
	}
	// Shoelace Formula, n var case
	double polyArea(ArrayList<Point> points) {
		double area = 0.0;
		points.add(points.get(0));
		for (int i=1; i<points.size(); ++i) {
			area += points.get(i-1).x*points.get(i).y-points.get(i-1).y*points.get(i).x;
		}
		area *= 0.5;
		points.remove(points.size()-1);
		return area;
	}
	// Heron's Formula
	double heronArea(Point A, Point B, Point C) {
		double a = B.dist(C), b = C.dist(A), c = A.dist(B);
		double s = (a+b+c)*0.5;
		return Math.sqrt(s*(s-a)*(s-b)*(s-c));
	}
	// determines if line _segments_ A1-A2 and B1-B2 intersect
	boolean segmentIntersect(Point A1, Point A2, Point B1, Point B2) {
		if (collinear(A1, A2, B1) || collinear(A1, A2, B2) ||
			collinear(B1, B2, A1) || collinear(B1, B2, A2))
			return false;
		return (Left(A1, A2, B1) ^ Left(A1, A2, B2)) && 
			   (Left(B1, B2, A1) ^ Left(B1, B2, A2));
	}
	// true iff lines A1-A2 and B1-B2 intersect - p. 33
	boolean lineIntersect(Point A1, Point A2, Point B1, Point B2) {
		if (segmentIntersect(A1, A2, B1, B2)) {
			return true;
		}
		if (between(A1, A2, B1) || between(A1, A2, B2) ||
			between(B1, B2, A1) || between(B1, B2, A2))
			return true;
		return false;
	}
	// true iff A,B,C collinear
	boolean collinear(Point A, Point B, Point C) {
		return triArea(A, B, C) < 1e-10;
	}
	// true iff C is to the left of line A->B
	boolean Left(Point A, Point B, Point C) {
		return triArea(A, B, C) > 0;
	}
	// true iff C is on AB between A and B
	boolean between(Point A, Point B, Point C) {
		if (!collinear(A, B, C)) {
			return false;
		}
		if (Math.abs(A.x-B.x) > 1e-10) {
			return ((A.x <= C.x) && (C.x <= B.x)) ||
				((A.x >= C.x) && (C.x >= B.x));
		} else {
			return ((A.y <= C.y) && (C.y <= B.y)) ||
				((A.y >= C.y) && (C.y >= B.y));
		}
	}
	double angle(Point A, Point B, Point C) {
		double ab = Math.sqrt((B.x-A.x)*(B.x-A.x)+(B.y-A.y)*(B.y-A.y));   
		double bc = Math.sqrt((B.x-C.x)*(B.x-C.x)+(B.y-C.y)*(B.y-C.y));
		double ac = Math.sqrt((C.x-A.x)*(C.x-A.x)+(C.y-A.y)*(C.y-A.y)); 
		return Math.acos((bc*bc+ab*ab-ac*ac)/(2*bc*ab));
	}
}
