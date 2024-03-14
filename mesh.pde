Mesh op_mesh;
final float lambda = 0.6;
final float lambda1 = 0.6307;
final float lambda2 = -0.6731;

class Mesh {
  //////////////////////////////////////////////////////////////
  ////////////////////// Vertex Class //////////////////////////
  //////////////////////////////////////////////////////////////
  class Vertex {
    PVector v;
    Face f; // one adjacent face
    Vertex(float _x, float _y, float _z) {
      v = new PVector(_x, _y, _z);
      f = null;
    }
    void setFace(Face _f) { if (f == null) f = _f; }
  }
  //////////////////////////////////////////////////////////////
  /////////////////////// Edge Class ///////////////////////////
  //////////////////////////////////////////////////////////////
  class Edge {
    Face f; // face that contains me
    int i; // starting index in f.vs
    SortedKey skey;
    PVector edgeMid;
    Edge(int _i, Face _f){ i = _i; f = _f;}
    SortedKey sortedKey() {
      if (skey == null) {
        int j = (i+1) % f.k;
        int vi = f.vs[i], vj = f.vs[j];
        skey = new SortedKey(vi, vj);
      }
      return skey;
    }
  }
  //////////////////////////////////////////////////////////////
  /////////////////////// Face Class ///////////////////////////
  //////////////////////////////////////////////////////////////
  class Face {
    int k; // how many sides
    int[] vs; // vertex indices
    PVector normal, fColor, midV;
    Key[] opp; // store opposite edge key, opp[i] is the opposite edge of <vs[i], vs[i+1]>
    Face(int _k, int[] _vs){ 
      k = _k; 
      vs = _vs;
      opp = new Key[k];
    }
  }
  //////////////////////////////////////////////////////////////
  ////////////////////// Mesh Members //////////////////////////
  //////////////////////////////////////////////////////////////
  ArrayList<Vertex> verts; // vertices
  ArrayList<Face> faces; // faces
  HashMap<Key, Edge> edges; // directed edge, <vi, vj> -> Edge
  ArrayList<PVector> vNormals; // vertex normal
  
  Mesh() {
    verts = new ArrayList<>();
    faces = new ArrayList<>();
    edges = new HashMap<>();
    vNormals = new ArrayList<>();
  }
  int nv() { return verts.size(); }
  int nf() { return faces.size(); }
  //////////////////////////////////////////////////////////////
  //////////////////////// Mesh Utils //////////////////////////
  //////////////////////////////////////////////////////////////
  void addVertex(float x, float y, float z) { 
    verts.add(new Vertex(x, y, z)); 
  }
  void addFace(int k, int[] v_indices) {
    Face f = new Face(k, v_indices);
    faces.add(f);
    for (int i = 0; i < k; i++) {
      // Set vertice's face
      verts.get(v_indices[i]).setFace(f);
      // Create Edge
      Edge e = new Edge(i, f);
      int j = (i + 1) % k;
      edges.put(new Key(v_indices[i], v_indices[j]), e);
      // Set Face opposite edge keys
      f.opp[i] = new Key(v_indices[j], v_indices[i]);
    }
    // face color
    if (f.fColor == null)
      f.fColor = new PVector(random(256), random(256), random(256));
  }
  // face mid point
  PVector facePoint(Face f) {
    if (f.midV != null) return f.midV;
    PVector m = new PVector();
    for (int i = 0; i < f.k; i++)
      m.add(verts.get(f.vs[i]).v);
    f.midV = m.div(f.k);
    return m;
  }
  // compute face normal
  PVector fNormal(Face f) {
    if (f.normal == null)
      f.normal = PVector.sub(verts.get(f.vs[1]).v, verts.get(f.vs[0]).v)
           .cross(PVector.sub(verts.get(f.vs[2]).v, verts.get(f.vs[0]).v))
           .normalize();
    return f.normal;
  }
  // compute vertex normal
  void computeVertexNormal() {
    for (int i = 0; i < nv(); i++) {
      PVector n = new PVector();
      // loop over v's adjacent faces
      Edge eStart = startFromVertex(i, verts.get(i));
      Edge e = eStart;
      do {
        n.add(fNormal(e.f));
        e = swing(e);
      } while (e != eStart);
      vNormals.add(n.normalize());
    }
  }
  // normalize vertex to sphere
  void normalize() {
    for (Vertex v: verts) v.v.normalize();
  }
  //////////////////////////////////////////////////////////////
  /////////////////////////// Dual /////////////////////////////
  //////////////////////////////////////////////////////////////
  Mesh dual() {
    Mesh dualMesh = new Mesh();
    HashMap<Face, Integer> f2vId = new HashMap();
    for (int i = 0; i < nv(); i++) {
      ArrayList<Integer> li = new ArrayList<Integer>(); //<>//
      Edge eStart = startFromVertex(i, verts.get(i));
      Edge e = eStart;
      do {
        if (!f2vId.containsKey(e.f)) {
          PVector v = facePoint(e.f);
          f2vId.put(e.f, dualMesh.nv());
          dualMesh.addVertex(v.x, v.y, v.z);
        }
        li.add(f2vId.get(e.f));
        e = unswing(e); // CCW
      } while (e != eStart);
      int[] faceInd = li.stream().mapToInt(Integer::valueOf).toArray();
      dualMesh.addFace(faceInd.length, faceInd);
    }
    dualMesh.computeVertexNormal();
    return dualMesh;
  }
  // move vertice along normal
  void noised() {
    for (int i = 0; i < nv(); i++) {
      Vertex vert = verts.get(i);
      float r = random(-0.1, 0.1);
      vert.v.x += vNormals.get(i).x * r;
      vert.v.y += vNormals.get(i).y * r;
      vert.v.z += vNormals.get(i).z * r;
    }
  }
  //////////////////////////////////////////////////////////////
  ///////////////////////// MidPoint ///////////////////////////
  //////////////////////////////////////////////////////////////
  Mesh subMid() {
    HashMap<SortedKey, Integer> edgeToVertex = new HashMap<>();
    
    Mesh subMesh = new Mesh();
    for (int i = 0; i < nv(); i++) {
      Vertex p = verts.get(i);
      Edge eStart = startFromVertex(i, p);
      Edge e = eStart;
      // 1. add myself
      int pIndex = subMesh.nv();
      subMesh.addVertex(p.v.x,p.v.y,p.v.z);
      do {
        int[] li = new int[3];
        li[0] = pIndex;
        // 2. edge middle
        PVector eMid1 = edgeMiddle(e);
        if (!edgeToVertex.containsKey(e.sortedKey())) {
          edgeToVertex.put(e.sortedKey(), subMesh.nv());
          subMesh.addVertex(eMid1.x,eMid1.y,eMid1.z);
        }
        li[1] = edgeToVertex.get(e.sortedKey());
        e = unswing(e);
        // 3. edge middle
        PVector eMid2 = edgeMiddle(e);
        if (!edgeToVertex.containsKey(e.sortedKey())) {
          edgeToVertex.put(e.sortedKey(), subMesh.nv());
          subMesh.addVertex(eMid2.x,eMid2.y,eMid2.z);
        }
        li[2] = edgeToVertex.get(e.sortedKey());
        // add face
        subMesh.addFace(li.length, li);
      } while (e != eStart);
    }
    // edge middle face
    for (int i = 0; i < nf(); i++) {
      Face f = faces.get(i);
      int[] li = new int[f.k];
      Edge eStart = edges.get(new Key(f.vs[0], f.vs[1]));
      Edge e = eStart;
      int cnt = 0;
      do {
        li[cnt++] = edgeToVertex.get(e.sortedKey());
        e = next(e);
      } while (e != eStart);
      subMesh.addFace(li.length, li);
    }
    // projection
    subMesh.normalize();
    // compute normal
    subMesh.computeVertexNormal();
    return subMesh;
  }
  PVector edgeMiddle(Edge e) {
    if (e.edgeMid != null) return e.edgeMid;
    int i = e.i, j = (i+1) % e.f.k;
    int vi = e.f.vs[i], vj = e.f.vs[j];
    Vertex Vi = verts.get(vi), Vj = verts.get(vj);
    PVector m = PVector.add(Vi.v, Vj.v).mult(0.5);
    e.edgeMid = m;
    return m;
  }
  //////////////////////////////////////////////////////////////
  ///////////////////////// Catmull ////////////////////////////
  //////////////////////////////////////////////////////////////
  Mesh subCat() {
    HashMap<SortedKey, Integer> edgeToVertex = new HashMap<>();
    HashMap<Face, Integer> faceToVertex = new HashMap<>();
    HashMap<SortedKey, PVector> edgePCache = new HashMap<>();

    Mesh subMesh = new Mesh();
    // Refine Topology
    for (int i = 0; i < nv(); i++) {
      Vertex p = verts.get(i);
      Edge eStart = startFromVertex(i, p);
      Edge e = eStart;
      PVector eSum = new PVector();
      PVector fSum = new PVector();
      int n = 0;
      // 1. add myself
      int pIndex = subMesh.nv();
      subMesh.addVertex(p.v.x,p.v.y,p.v.z);
      do {
        n++;
        // 1.1 add to new face
        int[] li = new int[4];
        li[0] = pIndex;
        // 2. edge point
        PVector edgeP1 = edgePoint(e, edgePCache);
        eSum.add(edgeMiddle(e)); // ! this is so fking annoying
        if (!edgeToVertex.containsKey(e.sortedKey())) {
          edgeToVertex.put(e.sortedKey(), subMesh.nv());
          subMesh.addVertex(edgeP1.x, edgeP1.y, edgeP1.z);
        }
        li[1] = edgeToVertex.get(e.sortedKey());
        // 3. face point
        PVector faceP = facePoint(e.f);
        fSum.add(faceP);
        if (!faceToVertex.containsKey(e.f)) {
          faceToVertex.put(e.f, subMesh.nv());
          subMesh.addVertex(faceP.x, faceP.y, faceP.z);
        }
        li[2] = faceToVertex.get(e.f);
        // ------------------------------------
        e = unswing(e); // CCW
        // ------------------------------------
        // 4. ccw edge point
        PVector edgeP2 = edgePoint(e, edgePCache);
        if (!edgeToVertex.containsKey(e.sortedKey())) {
          edgeToVertex.put(e.sortedKey(), subMesh.nv());
          subMesh.addVertex(edgeP2.x, edgeP2.y, edgeP2.z);
        }
        li[3] = edgeToVertex.get(e.sortedKey());
        // form face
        subMesh.addFace(li.length, li);
      } while (e != eStart);
      // Refine Geometry
      float coe = 1.0 / (float)n;
      eSum.mult(coe);
      fSum.mult(coe);
      PVector newP = new PVector();
      newP.add(PVector.mult(p.v ,(n-3) * coe));
      newP.add(eSum.mult(2 * coe));
      newP.add(fSum.mult(coe));
      subMesh.verts.get(pIndex).v = newP;
    }
    subMesh.computeVertexNormal();
    return subMesh;
  }
  PVector edgePoint(Edge e, HashMap<SortedKey, PVector> eCache) {
    if (eCache.containsKey(e.sortedKey())) return eCache.get(e.sortedKey());
    int i = e.i;
    int j = (i+1) % e.f.k;
    int vi = e.f.vs[i];
    int vj = e.f.vs[j];
    PVector Vi = verts.get(vi).v;
    PVector Vj = verts.get(vj).v;
    PVector fMid1 = facePoint(e.f);
    PVector fMid2 = facePoint(opp(e).f);
    PVector m = new PVector();
    m.add(Vi);
    m.add(Vj);
    m.add(fMid1);
    m.add(fMid2);
    m.div(4);
    eCache.put(e.sortedKey(), m);
    return m;
  }
  //////////////////////////////////////////////////////////////
  ///////////////////////// Smoothing //////////////////////////
  //////////////////////////////////////////////////////////////
  Mesh laplacian(int t) {
    PVector[] newVs1 = new PVector[nv()];
    PVector[] newVs2 = new PVector[nv()];
    for (int i = 0; i < nv(); i++) 
      newVs1[i] = verts.get(i).v;
    while (t-- > 0) {
      for (int i = 0; i < nv(); i++) {
        PVector oldP = newVs1[i];
        Edge eStart = startFromVertex(i, verts.get(i));
        Edge e = eStart;
        PVector newP = new PVector();
        int n = 0;
        do {
          n++;
          int j = (e.i + 1) % e.f.k;
          newP.add(newVs1[e.f.vs[j]]);
          e = swing(e);
        } while (e != eStart);
        newP.div(n).sub(oldP).mult(lambda).add(oldP);
        newVs2[i] = newP;
      }
      PVector[] tmp = newVs1;
      newVs1 = newVs2;
      newVs2 = tmp;
    }
    Mesh newMesh = new Mesh();
    for (int i = 0; i < nv(); i++) 
      newMesh.addVertex(newVs1[i].x, newVs1[i].y, newVs1[i].z);
    for (Face f: faces) {
      int[] ind = f.vs.clone();
      newMesh.addFace(f.k, ind);
    }
    newMesh.computeVertexNormal();
    return newMesh;
  }
  Mesh taubin(int t) {
    PVector[] newVs1 = new PVector[nv()];
    PVector[] newVs2 = new PVector[nv()];
    for (int i = 0; i < nv(); i++) 
      newVs1[i] = verts.get(i).v;
    while (t-- > 0) {
      for (int i = 0; i < nv(); i++) {
        PVector oldP = newVs1[i];
        Edge eStart = startFromVertex(i, verts.get(i));
        Edge e = eStart;
        PVector newP = new PVector();
        int n = 0;
        do {
          n++;
          int j = (e.i + 1) % e.f.k;
          newP.add(newVs1[e.f.vs[j]]);
          e = swing(e);
        } while (e != eStart);
        newP.div(n).sub(oldP).mult(lambda1).add(oldP);
        newVs2[i] = newP;
      }
      for (int i = 0; i < nv(); i++) {
        PVector oldP = newVs2[i];
        Edge eStart = startFromVertex(i, verts.get(i));
        Edge e = eStart;
        PVector newP = new PVector();
        int n = 0;
        do {
          n++;
          int j = (e.i + 1) % e.f.k;
          newP.add(newVs2[e.f.vs[j]]);
          e = swing(e);
        } while (e != eStart);
        newP.div(n).sub(oldP).mult(lambda2).add(oldP);
        newVs1[i] = newP;
      }
    }
    Mesh newMesh = new Mesh();
    for (int i = 0; i < nv(); i++) 
      newMesh.addVertex(newVs1[i].x, newVs1[i].y, newVs1[i].z);
    for (Face f: faces) {
      int[] ind = f.vs.clone();
      newMesh.addFace(f.k, ind);
    }
    newMesh.computeVertexNormal();
    return newMesh;
  }
  //////////////////////////////////////////////////////////////
  //////////////////////// Edge Utils //////////////////////////
  //////////////////////////////////////////////////////////////
  Edge next(Edge e) {
    int i = (e.i+1) % e.f.k;
    int j = (i+1) % e.f.k;
    int vi = e.f.vs[i];
    int vj = e.f.vs[j];
    return edges.get(new Key(vi, vj));
  }
  Edge prev(Edge e) {
    int i = (e.i-1+e.f.k) % e.f.k;
    int j = e.i;
    int vi = e.f.vs[i];
    int vj = e.f.vs[j];
    return edges.get(new Key(vi, vj));
  }
  Edge opp(Edge e) {
    return edges.get(e.f.opp[e.i]);
  }
  Edge swing(Edge e) {
    return next(opp(e));
  }
  Edge unswing(Edge e) {
    return opp(prev(e));
  }
  Edge startFromVertex(int vi, Vertex v) {
    Face f = v.f;
    int j;
    for (j = 0; j < f.k; j++)
      if (f.vs[j] == vi)
        break;
    int vj = f.vs[(j + 1) % f.k];
    
    return edges.get(new Key(vi, vj));
  }
  //////////////////////////////////////////////////////////////
  /////////////////////// Render Utils /////////////////////////
  //////////////////////////////////////////////////////////////
  void render(boolean randomColor) {
    for (Face f: faces) {
      beginShape();
      if (randomColor) fill(f.fColor.x, f.fColor.y, f.fColor.z);
      for (int i: f.vs) {
        if (use_v_normal) normal(vNormals.get(i).x, vNormals.get(i).y, vNormals.get(i).z);
        vertex(verts.get(i).v.x, verts.get(i).v.y, verts.get(i).v.z);
      }
      endShape(CLOSE);
    }
  }
  Edge op_edge;
  void renderEdge(boolean flag) {
    if (!flag) return;
    if (op_edge == null) op_edge = (Edge)edges.values().toArray()[0];
    Edge e = op_edge;
    // get mid point
    int i = e.i;
    int j = (i+1) % e.f.k;
    int vi = e.f.vs[i];
    int vj = e.f.vs[j];
    Vertex Vi = verts.get(vi), Vj = verts.get(vj);
    PVector dir = PVector.sub(Vj.v, Vi.v);
    PVector w = fNormal(e.f).cross(dir).normalize();
    PVector mid = PVector.add(Vi.v, PVector.mult(dir, 0.50));
    PVector p1 = PVector.add(mid, PVector.mult(w, 0.04));
    dir.normalize();
    PVector p0 = PVector.add(p1, PVector.mult(dir, -0.04));
    PVector p2 = PVector.add(p1, PVector.mult(dir, 0.04));
    noStroke();
    fill(204, 102, 0);
    pushMatrix();
    translate(p0.x, p0.y, p0.z);
    sphere(0.04);
    popMatrix();
    pushMatrix();
    translate(p1.x, p1.y, p1.z);
    sphere(0.03);
    popMatrix();
    pushMatrix();
    translate(p2.x, p2.y, p2.z);
    sphere(0.025);
    popMatrix();
  }
  void printMe() {
    println("vertex " + nv());
    println("face " + nf());
    for (Vertex v: verts) {
      print(v.v.x + " ");
      print(v.v.y + " ");
      print(v.v.z);
      println();
    }
    for (Face f: faces) {
      print(f.k);
      for (int i = 0; i < f.k; i++)
        print(" " + f.vs[i]);
      println();
    }
  }
}
