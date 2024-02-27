#include "student_code.h"
#include "mutablePriorityQueue.h"
#include <iostream>

#define nl "\n"

using namespace std;

namespace CGL
{


	Vector2D lerp2D(Vector2D a, Vector2D b, double d) {

		return d * a + (1 - d) * b;

	}

	Vector3D lerp3D(Vector3D a, Vector3D b, double d) {

		return d * a + (1 - d) * b;

	}


	/**
	 * Evaluates one step of the de Casteljau's algorithm using the given points and
	 * the scalar parameter t (class member).
	 *
	 * @param points A vector of points in 2D
	 * @return A vector containing intermediate points or the final interpolated vector
	 */
	std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const& points)
	{
		// TODO Part 1.

		vector<Vector2D> vec(points.size() - 1);

		for (int i = 0; i < vec.size(); i++) {
			Vector2D cur = points.at(i);
			Vector2D next = points.at(i + 1);

			vec.at(i) = lerp2D(cur, next, t);

		}

		return vec;
		//return std::vector<Vector2D>();
	}

	/**
	 * Evaluates one step of the de Casteljau's algorithm using the given points and
	 * the scalar parameter t (function parameter).
	 *
	 * @param points    A vector of points in 3D
	 * @param t         Scalar interpolation parameter
	 * @return A vector containing intermediate points or the final interpolated vector
	 */
	std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const& points, double t) const
	{
		// TODO Part 2.

		vector<Vector3D> vec(points.size() - 1);

		for (int i = 0; i < vec.size(); i++) {
			Vector3D cur = points.at(i);
			Vector3D next = points.at(i + 1);

			vec.at(i) = lerp3D(cur, next, t);

		}

		return vec;
	}

	/**
	 * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
	 *
	 * @param points    A vector of points in 3D
	 * @param t         Scalar interpolation parameter
	 * @return Final interpolated vector
	 */
	Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const& points, double t) const
	{
		// TODO Part 2.

		vector<Vector3D> vec(points.begin(), points.end());

		while (vec.size() > 1) {
			vec = evaluateStep(vec, t);
		}

		return vec.at(0);

		return Vector3D();
	}

	/**
	 * Evaluates the Bezier patch at parameter (u, v)
	 *
	 * @param u         Scalar interpolation parameter
	 * @param v         Scalar interpolation parameter (along the other axis)
	 * @return Final interpolated vector
	 */
	Vector3D BezierPatch::evaluate(double u, double v) const
	{
		int n = controlPoints.size();

		vector<Vector3D> vec(n);

		for (int i = 0; i < n; i++) {
			vec.at(i) = evaluate1D(controlPoints.at(i), u);
		}

		return evaluate1D(vec, v);


		// TODO Part 2
	}

	Vector3D Vertex::normal(void) const
	{
		// TODO Part 3.
		// Returns an approximate unit normal at this vertex, computed by
		// taking the area-weighted average of the normals of neighboring
		// triangles, then normalizing.
		  /*cout << position << " ";*/
		HalfedgeIter h = _halfedge;

		Vector3D total = Vector3D();

		do
		{
			// don't count boundary loops
			if (!h->face()->isBoundary())
			{
				/*cout << h->face()->normal() << " ";*/
				total += h->face()->normal();
			}

			// move to the next halfedge around the vertex
			h = h->twin()->next();
		} while (h != _halfedge); // done iterating over halfedges

		total.normalize();

		//cout << total;
		return total;

	}

	EdgeIter HalfedgeMesh::flipEdge(EdgeIter e0)
	{

		Edge* e = e0->getEdge();

		HalfedgeIter bc = e->halfedge();
		HalfedgeIter cb = bc->getHalfedge()->twin();

		Halfedge* h = bc->getHalfedge();

		Halfedge* twin = cb->getHalfedge();

		HalfedgeIter ca = h->next();
		HalfedgeIter ab = ca->next();

		HalfedgeIter bd = twin->next();
		HalfedgeIter dc = bd->next();

		FaceIter abc = h->face();
		FaceIter cbd = twin->face();

		VertexIter b = h->vertex();
		VertexIter c = twin->vertex();
		VertexIter a = h->next()->getHalfedge()->next()->getHalfedge()->vertex();
		VertexIter d = twin->next()->getHalfedge()->next()->getHalfedge()->vertex();



		if (abc->getFace()->isBoundary() || cbd->getFace()->isBoundary() || e0->getEdge()->isBoundary()) {
			return e0;
		}


		ab->setNeighbors(bd, ab->getHalfedge()->twin(), a, ab->edge(), abc);
		bd->setNeighbors(bc, bd->getHalfedge()->twin(), b, bd->edge(), abc);
		bc->setNeighbors(ab, bc->getHalfedge()->twin(), d, bc->edge(), abc);

		dc->setNeighbors(ca, dc->getHalfedge()->twin(), d, dc->edge(), cbd);
		ca->setNeighbors(cb, ca->getHalfedge()->twin(), c, ca->edge(), cbd);
		cb->setNeighbors(dc, cb->getHalfedge()->twin(), a, cb->edge(), cbd);


		abc->getFace()->halfedge() = bc;
		a->getVertex()->halfedge() = cb;
		b->getVertex()->halfedge() = bd;
		c->getVertex()->halfedge() = ca;
		d->getVertex()->halfedge() = dc;

		abc->getFace()->halfedge() = bc;
		cbd->getFace()->halfedge() = cb;



		return cb->edge();


		// TODO Part 4.
		// This method should flip the given edge and return an iterator to the flipped edge.
		return EdgeIter();
	}



	VertexIter HalfedgeMesh::splitEdge(EdgeIter e0)
	{
		HalfedgeIter bc = e0->getEdge()->halfedge();
		HalfedgeIter cb = bc->getHalfedge()->twin();

		FaceIter abc = bc->getHalfedge()->face();
		FaceIter cbd = cb->getHalfedge()->face();

		if (abc->getFace()->isBoundary() || cbd->getFace()->isBoundary() || e0->getEdge()->isBoundary()) {

			// inits
			HalfedgeIter ca = bc->getHalfedge()->next();
			HalfedgeIter ab = ca->next();

			VertexIter b = bc->getHalfedge()->vertex();
			VertexIter c = cb->getHalfedge()->vertex();
			VertexIter a = bc->getHalfedge()->next()->getHalfedge()->next()->getHalfedge()->vertex();

			VertexIter m = newVertex();
			m->getVertex()->position = (b->getVertex()->position + c->getVertex()->position) * 0.5;

			FaceIter amc = newFace();
			FaceIter bma = newFace();

			EdgeIter am_e = newEdge();
			EdgeIter bm_e = newEdge();
			EdgeIter cm_e = newEdge();


			HalfedgeIter am = newHalfedge(), ma = newHalfedge();
			HalfedgeIter bm = newHalfedge(), mb = newHalfedge();
			HalfedgeIter mc = newHalfedge();

			m->getVertex()->halfedge() = mc;
			b->getVertex()->halfedge() = bm;
			c->getVertex()->halfedge() = ca;

			amc->getFace()->halfedge() = am;
			bma->getFace()->halfedge() = bm;

			am_e->getEdge()->halfedge() = am;
			bm_e->getEdge()->halfedge() = bm;
			cm_e->getEdge()->halfedge() = mc;

			//
			ma->setNeighbors(ab, am, m, am_e, bma);
			mb->setNeighbors(cb->next(), bm, m, bm_e, cbd);
			mc->setNeighbors(ca, cb, m, cm_e, amc);


			am->setNeighbors(mc, ma, a, am_e, amc);
			bm->setNeighbors(ma, mb, b, bm_e, bma);
			cb->setNeighbors(mb, mc, c, cm_e, cbd);

			ab->setNeighbors(bm, ab->getHalfedge()->twin(), a, ab->getHalfedge()->edge(), bma);
			ca->setNeighbors(am, ca->getHalfedge()->twin(), c, ca->getHalfedge()->edge(), amc);
			//make cb the new cm.

			deleteEdge(e0);
			deleteHalfedge(bc);
			deleteFace(abc);

			return VertexIter();
		}


		HalfedgeIter ca = bc->getHalfedge()->next();
		HalfedgeIter ab = ca->next();
		HalfedgeIter bd = cb->getHalfedge()->next();
		HalfedgeIter dc = bd->next();

		VertexIter b = bc->getHalfedge()->vertex();
		VertexIter c = cb->getHalfedge()->vertex();
		VertexIter a = bc->getHalfedge()->next()->getHalfedge()->next()->getHalfedge()->vertex();
		VertexIter d = cb->getHalfedge()->next()->getHalfedge()->next()->getHalfedge()->vertex();

		VertexIter m = newVertex();
		m->getVertex()->position = (b->getVertex()->position + c->getVertex()->position) * 0.5;
		m->getVertex()->newPosition = e0->getEdge()->newPosition;
		
		m->getVertex()->isNew = true;


		FaceIter amc = newFace();
		FaceIter cmd = newFace();
		FaceIter dmb = newFace();
		FaceIter bma = newFace();


		EdgeIter am_e = newEdge();
		EdgeIter bm_e = newEdge();
		//EdgeIter cm_e = newEdge();
		EdgeIter dm_e = newEdge();

		am_e->getEdge()->isNew = true;
		dm_e->getEdge()->isNew = true;
		bm_e->getEdge()->isNew = false;
		e0->getEdge()->isNew = false;

		am_e->getEdge()->alrSplit = true;
		dm_e->getEdge()->alrSplit = true;
		bm_e->getEdge()->alrSplit = true;
		e0->getEdge()->alrSplit = true;

		HalfedgeIter am = newHalfedge(), ma = newHalfedge();
		HalfedgeIter bm = newHalfedge(), mb = newHalfedge();
		HalfedgeIter cm = newHalfedge(), mc = newHalfedge();
		HalfedgeIter dm = newHalfedge(), md = newHalfedge();

		m->getVertex()->halfedge() = mc;
		b->getVertex()->halfedge() = bm;
		c->getVertex()->halfedge() = cm;

		amc->getFace()->halfedge() = am;
		cmd->getFace()->halfedge() = cm;
		dmb->getFace()->halfedge() = dm;
		bma->getFace()->halfedge() = bm;

		am_e->getEdge()->halfedge() = am;
		bm_e->getEdge()->halfedge() = bm;
		e0->getEdge()->halfedge() = cm;
		dm_e->getEdge()->halfedge() = dm;

		am->setNeighbors(mc, ma, a, am_e, amc);
		bm->setNeighbors(ma, mb, b, bm_e, bma);
		cm->setNeighbors(md, mc, c, e0, cmd);
		dm->setNeighbors(mb, md, d, dm_e, dmb);


		ma->setNeighbors(ab, am, m, am_e, bma);
		mb->setNeighbors(bd, bm, m, bm_e, dmb);
		mc->setNeighbors(ca, cm, m, e0, amc);
		md->setNeighbors(dc, dm, m, dm_e, cmd);

		ab->setNeighbors(bm, ab->getHalfedge()->twin(), a, ab->getHalfedge()->edge(), bma);
		bd->setNeighbors(dm, bd->getHalfedge()->twin(), b, bd->getHalfedge()->edge(), dmb);
		dc->setNeighbors(cm, dc->getHalfedge()->twin(), d, dc->getHalfedge()->edge(), cmd);
		ca->setNeighbors(am, ca->getHalfedge()->twin(), c, ca->getHalfedge()->edge(), amc);

		
		//deleteEdge(e0);
		deleteHalfedge(bc); deleteHalfedge(cb);
		deleteFace(abc); deleteFace(cbd);
		return m;
	}



	void MeshResampler::upsample(HalfedgeMesh& mesh)
	{
		// TODO Part 6.
		// This routine should increase the number of triangles in the mesh using Loop subdivision.
		// One possible solution is to break up the method as listed below.

		// 1. Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
		// and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
		// a vertex of the original mesh.
		cout << "here\n";
		for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
			double n = (double)v->getVertex()->degree();
			double u = 3.0 / (8.0 * n);

			if (n == 3.0) {
				u = 3.0 / 16.0;

				
			}


			v->getVertex()->newPosition = (1.0 - n * u) * v->getVertex()->position + (u * v->getVertex()->neighborSum());

			v->getVertex()->isNew = false;

		}

		




			// 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.
		for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {

			HalfedgeIter h = e->getEdge()->halfedge();
			HalfedgeIter twin = h->getHalfedge()->twin();

			VertexIter a = h->getHalfedge()->next()->getHalfedge()->next()->getHalfedge()->vertex();
			VertexIter d = twin->getHalfedge()->next()->getHalfedge()->next()->getHalfedge()->vertex();
			VertexIter b = h->getHalfedge()->vertex();
			VertexIter c = twin->getHalfedge()->vertex();

			e->getEdge()->newPosition = (a->getVertex()->position + d->getVertex()->position) * 1.0 / 8.0 + (b->getVertex()->position + c->getVertex()->position) * 3.0 / 8.0;
			//cout << e->getEdge() << a->getVertex()->position << " " << b->getVertex()->position << " " << c->getVertex()->position << " " << d->getVertex()->position << nl;
			e->getEdge()->alrSplit = false;
			e->getEdge()->isNew = false;
		}
			// 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
			// information about which subdivide edges come from splitting an edge in the original mesh, and which edges
			// are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
			// the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)
		
		EdgeIter end = mesh.edgesEnd();
		int ee = 0;

		cout << mesh.nEdges() << nl;

		for (EdgeIter e = mesh.edgesBegin(); e != end; e++) {

			if (!e->getEdge()->alrSplit) {
				mesh.splitEdge(e);
			}

		}
			// 4. Flip any new edge that connects an old and new vertex.
		
		cout << mesh.nEdges() << nl;

		for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {

			HalfedgeIter ab = e->getEdge()->halfedge();
			HalfedgeIter ba = ab->getHalfedge()->twin();

			VertexIter a = ab->getHalfedge()->vertex();
			VertexIter b = ba->getHalfedge()->vertex();

			if ((a->getVertex()->isNew && !b->getVertex()->isNew) || (b->getVertex()->isNew && !a->getVertex()->isNew)) {
				if (e->getEdge()->isNew) {
					mesh.flipEdge(e);
				}
				
			}


		}



			// 5. Copy the new vertex positions into final Vertex::position.

		for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
			v->getVertex()->position = v->getVertex()->newPosition;
		}



		}
	}
