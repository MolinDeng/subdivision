// create a mesh
void read_mesh(String filename)
{
  String[] words;
  
  String lines[] = loadStrings(filename);
  
  words = split (lines[0], " ");
  int num_vertices = int(words[1]);
  
  words = split (lines[1], " ");
  int num_faces = int(words[1]);
  
  Mesh mesh = new Mesh();
  
  // read in the vertices
  for (int i = 0; i < num_vertices; i++) {
    words = split (lines[i+2], " ");
    float x = float(words[0]);
    float y = float(words[1]);
    float z = float(words[2]);
    mesh.addVertex(x, y, z);
  }
  
  // read in the faces
  for (int i = num_faces-1; i >= 0; i--) { // CCW
  //for (int i = 0; i < num_faces; i++) { // CW
    int j = i + num_vertices + 2;
    words = split (lines[j], " ");
    // get the number of vertices for this face
    int nverts = int(words[0]);
    // get all of the vertex indices
    int[] vertices = new int[nverts];
    for (int k = 1; k <= nverts; k++) 
      vertices[k-1] = int(words[k]);
    mesh.addFace(nverts, vertices);
  }
  op_mesh = mesh;
  op_mesh.computeVertexNormal();
}
