//// #include "./utils/WriteBMP.h"
//#include "./utils/MyQuaternion.h"
//
//#include <glm/glm.hpp>
//#include <glm/gtc/matrix_transform.hpp>
//#include <glm/gtc/type_ptr.hpp>
//#include <igl/readOFF.h>
//
//#include <igl/lexicographic_triangulation.h>
//#include <igl/opengl/glfw/Viewer.h>
//#include <igl/triangle/triangulate.h>
//#include <igl/gaussian_curvature.h>
//#include <igl/massmatrix.h>
//#include <igl/invert_diag.h>
//#include <igl/read_triangle_mesh.h>
//#include <igl/png/writePNG.h>
//#include <igl/PI.h>
//#include <Eigen/Geometry>
//// embree
//#include <igl/embree/EmbreeRenderer.h>
//
//#include <iostream>
//
///* mesh input */
//Eigen::MatrixXd inputV1;
//Eigen::MatrixXi inputF1;
//Eigen::MatrixXd inputC1;
//const char *mesh_pathname;
//std::string storage_pathname;
//
//// Default Values of Variables
//double normal_x = 0.0;
//double normal_y = 1.0;
//double normal_z = 0.0;
//Eigen::Vector3d old_normal(1.0,0.0,0.0);
//Eigen::Vector3d new_normal(0.0,1.0,0.0);
//Eigen::Vector3d rot_axis;
//float rot_angle = 0.0f;
//int num_slices = 0;
//double dis_slices = 100.0;
//double plane_xcord = 0; // default transations
//double plane_ycord = 0;
//double plane_zcord = 0;
//double plane_offset = 0;
//bool contour = false;
//bool image = false;
//bool color = false;
//
///* Cutting Plane coorinate */
//double max_value, min_value;
//// double plane[] = { // default plane
////     -0.01f, -1.0f, 1.0f,
////     -0.01f, -1.0f, -1.0f,
////     -0.01f, 1.0f, 1.0f,
////     -0.01f, 1.0f, -1.0f};
//double plane[] = { // default plane
//    -1.0f,  1.0f, -0.01f,
//    -1.0f, -1.0f, -0.01f,
//     1.0f,  1.0f, -0.01f,
//     1.0f, -1.0f, -0.01f};
//double temp_plane[12]; // creating a copy for modification
//unsigned int planeind[] = {
//    0, 1, 2,
//    1, 3, 2};
//
//// Default drawing bool values
//bool draw_mesh = true; // start with true;
//bool draw_wire = false;
//bool draw_plane = false;
//bool draw_contour = false;
//bool allow_cursor_movement = false;
//
///* for input mesh */
//struct Vertex
//{
//    glm::vec3 Position;    // Position
//    glm::vec3 vertexColor; // Vertex color
//};
//
///* contour */
//std::vector<Eigen::RowVector3d> outputVertices;
//std::vector<Eigen::RowVector2d> outputEdges;
//std::vector<Eigen::RowVector3d> outputColor;
//std::vector<Vertex> slices;
//
//// function declarations
//bool point_onPlane(Eigen::Vector3d &vec1, Eigen::Vector3d &nomPL, Eigen::Vector3d &ptPL);
//bool plane_triangleIntersect(Eigen::Vector3d &vec1, Eigen::Vector3d &vec2, Eigen::Vector3d &vec3, Eigen::Vector3d &nomPL, Eigen::Vector3d &ptPL);
//std::vector<Eigen::Vector3d> plane_edgeIntersect(Eigen::Vector3d &vec1, Eigen::Vector3d &vec2, Eigen::Vector3d &vec3, Eigen::Vector3d &nomPL, Eigen::Vector3d &ptPL);
//void compute_intersection(Eigen::MatrixXd &mesh_ver, Eigen::MatrixXi &mesh_fac, Eigen::MatrixXd &mesh_col, Eigen::Vector3d &normal_PL, Eigen::Vector3d &point_PL);
//void writePLY(std::vector<Eigen::RowVector3d> VERTICES, std::vector<Eigen::RowVector3d> COLOR, std::vector<Eigen::RowVector2d> EDGES, std::string directory);
//void writePNG_FILE(std::vector<Eigen::RowVector3d> VERTICES, std::vector<Eigen::RowVector3d> COLORS, std::vector<Eigen::RowVector2d> EDGES, std::string PNG_FILE);
//std::vector<Eigen::Vector3d> rayBB_intersection(Eigen::Vector3d &min_BB, Eigen::Vector3d &max_BB, Eigen::Vector3d &ray_direction);
//void updatePoints(std::vector<Eigen::RowVector3d> &outputVertices, Eigen::Vector4d &size);
//
//int main(int argc, char **argv)
//{
//
//    if (argc != 3)
//    {
//        std::cout << "Not enough arguments" << std::endl;
//        std::cout << "Usage: mesh_slicer [mesh_file_name.off] [storage_directory_location]" << std::endl;
//        return 1;
//    }
//    else
//    {
//
//        mesh_pathname = argv[1];
//        storage_pathname = argv[2];
//    }
//
//    /* ------------------------------Mesh Data to Load------------------------------------------------------ */
//    igl::readOFF(mesh_pathname, inputV1, inputF1, inputC1);
//
//    std::cout << " INPUT =====================================" << std::endl;
//    std::cout << " mesh_dir = " << mesh_pathname << std::endl;
//    std::cout << " storage_dir = " << storage_pathname << std::endl;
//    std::cout << " #Vertices in input_mesh = " << inputV1.rows() << std::endl;
//    std::cout << " #Faces in input_mesh = " << inputF1.rows() << std::endl;
//
//    if (inputC1.rows() > 0)
//    {
//        std::cout << " Vertices of input_mesh have color " << std::endl;
//        color = true;
//    }
//    else
//    {
//        std::cout << " Vertices of input_mesh have NO color " << std::endl;
//    }
//    std::cout << " =====================================" << std::endl;
//
//    std::vector<Vertex> vertices;
//    std::vector<size_t> indices;
//
//    min_value = inputV1.minCoeff(); // compute the overall min as we dont want to uniformaly scale the mesh down
//    max_value = inputV1.maxCoeff();
//    std::cout << "min_value : " << min_value << std::endl;
//    std::cout << "max_value : " << max_value << std::endl;
//
//    /* Walk through each of the mesh's vertices */
//    for (size_t i = 0; i < inputV1.rows(); i++)
//    {
//        Vertex vertex;
//        glm::vec3 vec;
//
//        // Positions --- Normalize each vertex
//        vec.x = 2 * ((inputV1(i, 0) - min_value) / (max_value - min_value)) - 1;
//        vec.y = 2 * ((inputV1(i, 1) - min_value) / (max_value - min_value)) - 1;
//        vec.z = 2 * ((inputV1(i, 2) - min_value) / (max_value - min_value)) - 1;
//        vertex.Position = vec;
//
//        if (color == true)
//        {
//            // VertexColor
//            vec.x = inputC1(i, 0);
//            vec.y = inputC1(i, 1);
//            vec.z = inputC1(i, 2);
//            vertex.vertexColor = vec;
//        }
//        else
//        {
//
//            vertex.vertexColor = glm::vec3(1.0f, 0.0f, 1.0f); // magenta color mesh
//        }
//        vertices.push_back(vertex);
//    }
//    // Walk through each of the mesh's faces
//    for (size_t i = 0; i < inputF1.rows(); i++)
//    {
//        indices.push_back(inputF1(i, 0));
//        indices.push_back(inputF1(i, 1));
//        indices.push_back(inputF1(i, 2));
//    }
//    // initialize temp_plane
//    for (int i = 0; i < 12; ++i)
//    {
//        temp_plane[i] = plane[i];
//    }
//
//    /* computing the angle of rotation for plane */
//    new_normal << normal_x, normal_y, normal_z;
//    new_normal.normalize(); // 0 0 1
//    old_normal.normalize(); // 1 0 0
//    if (old_normal.dot(new_normal) <= 1 && old_normal.dot(new_normal) >= -1 && old_normal != new_normal)
//    {
//        rot_angle = acos(old_normal.dot(new_normal));
//        rot_axis = old_normal.cross(new_normal);
//    }
//    else if (old_normal == new_normal)
//    {
//        rot_angle = 0;
//        rot_axis = old_normal;
//    }
//
//    contour = true;
//    double intersectDistance, inter_sliceDist;
//    int number_slices;
//    /* check which two pts does the ray (pt with new normal) intersect with the BB */
//    /* can be done using the min and max of the Bounding box only */
//    Eigen::Vector3d BB_min, BB_max;
//    BB_min = inputV1.colwise().minCoeff();
//    BB_max = inputV1.colwise().maxCoeff();
//
//    BB_min = BB_min.array() + 5; // slightly shrinking the bounding box to allow intersection
//    BB_max = BB_max.array() - 5;
//
//    std::vector<Eigen::Vector3d> intersectPts;
//    Eigen::Vector3d new_point;
//
//    intersectPts = rayBB_intersection(BB_min, BB_max, new_normal);
//
//    intersectDistance = sqrt(pow(intersectPts.at(0)(0) - intersectPts.at(1)(0), 2) + pow(intersectPts.at(0)(1) - intersectPts.at(1)(1), 2) + pow(intersectPts.at(0)(2) - intersectPts.at(1)(2), 2)); // distance
//
//    if (num_slices == 0)
//    {
//        inter_sliceDist = dis_slices / 4; // the ratio 100um -> 25um
//        number_slices = intersectDistance / inter_sliceDist;
//    }
//    else
//    {
//        number_slices = num_slices;
//        inter_sliceDist = intersectDistance / (number_slices - 1);
//    }
//
//    if (contour)
//    {
//        std::cout << "Generating contours----" << std::endl;
//    }
//    else
//    {
//        std::cout << "Generating images-----" << std::endl;
//    }
//
//    for (int i = 0; i < number_slices; i++)
//    { // generating the contours or images
//
//        new_point = intersectPts.at(0) + (i * inter_sliceDist) * new_normal;
//        compute_intersection(inputV1, inputF1, inputC1, new_normal, new_point); // slice the mesh
//
//        if (outputVertices.size() > 0)
//        { // there might be no vertices at all from the intersection
//
//            std::cout << "contour # " << i + 1 << std::endl;
//            if (contour)
//            {
//
//                //writePLY(outputVertices, outputColor, outputEdges, storage_pathname + "/" + std::to_string(i+1)+".ply");
//                 writePNG_FILE(outputVertices, outputColor, outputEdges, storage_pathname + "/" + std::to_string(i + 1) + ".png");
//            }
//            else
//            {
//
//                // Eigen::Vector4d sz;
//                // updatePoints(outputVertices, sz); // aligning the vertices to z axis for writing to bmp
//                // WriteBMP(ceil(sz(0)) + 50, ceil(sz(1)) + 50, outputVertices, outputColor, outputEdges, storage_pathname + "/" + std::to_string(i + 1) + ".bmp");
//            }
//        }
//
//        outputVertices.clear(); // clearing the old values
//        outputColor.clear();
//        outputEdges.clear();
//    }
//    std::cout << "Done----" << std::endl;
//
//    return 0;
//}
//
//void compute_intersection(Eigen::MatrixXd &mesh_ver, Eigen::MatrixXi &mesh_fac, Eigen::MatrixXd &mesh_col, Eigen::Vector3d &normal_PL, Eigen::Vector3d &point_PL)
//{
//
//    std::vector<Eigen::Vector3d> vertices_temp;
//
//    int j = 0; // counter for vertices and color
//    int k = 0; // counter for the edges
//
//    /* check whether triangles intersect the plane or not */
//    for (int i = 0; i < mesh_fac.rows(); i++)
//    {
//
//        Eigen::Vector3d vec1 = mesh_ver.row(mesh_fac(i, 0)); // vertices of triangle faces
//        Eigen::Vector3d vec2 = mesh_ver.row(mesh_fac(i, 1));
//        Eigen::Vector3d vec3 = mesh_ver.row(mesh_fac(i, 2));
//
//        if (plane_triangleIntersect(vec1, vec2, vec3, normal_PL, point_PL))
//        {
//            /* there is an intersection */
//
//            vertices_temp = plane_edgeIntersect(vec1, vec2, vec3, normal_PL, point_PL);
//            /* vertices_temp.rows() ==2 or 3 cant be anything else*/
//
//            if (vertices_temp.size() == 3)
//            {
//
//                outputVertices.push_back(vertices_temp.at(0));
//                outputVertices.push_back(vertices_temp.at(1));
//                outputVertices.push_back(vertices_temp.at(2));
//
//                Eigen::Vector2d t(j, j + 1);
//                outputEdges.push_back(t);
//                t << j + 1, j + 2;
//                outputEdges.push_back(t);
//                t << j + 2, j;
//                outputEdges.push_back(t);
//
//                outputColor.push_back(mesh_col.row(mesh_fac(i, 0)));
//                outputColor.push_back(mesh_col.row(mesh_fac(i, 0)));
//                outputColor.push_back(mesh_col.row(mesh_fac(i, 0)));
//
//                j = j + 3;
//                k = k + 3;
//            }
//            else if (vertices_temp.size() == 2)
//            {
//
//                outputVertices.push_back(vertices_temp.at(0));
//                outputVertices.push_back(vertices_temp.at(1));
//
//                Eigen::Vector2d t(j, j + 1);
//                outputEdges.push_back(t);
//
//                outputColor.push_back(mesh_col.row(mesh_fac(i, 0)));
//                outputColor.push_back(mesh_col.row(mesh_fac(i, 0)));
//
//                j = j + 2;
//                k = k + 1;
//            }
//        }
//    }
//}
//
//std::vector<Eigen::Vector3d> plane_edgeIntersect(Eigen::Vector3d &vec1, Eigen::Vector3d &vec2, Eigen::Vector3d &vec3, Eigen::Vector3d &nomPL, Eigen::Vector3d &ptPL)
//{
//
//    /* compute the intersection of line segments (traiangle edges) with plane */
//    std::vector<Eigen::Vector3d> intersects;
//
//    /* all three pts lie on plane */
//    if (point_onPlane(vec1, nomPL, ptPL) && point_onPlane(vec2, nomPL, ptPL) && point_onPlane(vec3, nomPL, ptPL))
//    {
//
//        intersects.push_back(vec1);
//        intersects.push_back(vec2);
//        intersects.push_back(vec3);
//        return intersects;
//    }
//
//    Eigen::Matrix<double, 4, 3> vertices; /* new matrix for looping those the edges */
//    vertices.row(0) = vec1;
//    vertices.row(1) = vec2;
//    vertices.row(2) = vec3;
//    vertices.row(3) = vec1;
//
//    for (int i = 0; i < 3; i++)
//    {
//
//        Eigen::Vector3d v1 = vertices.row(i);
//        Eigen::Vector3d v2 = vertices.row(i + 1);
//
//        if (nomPL.dot(v2 - v1) == 0)
//        {
//            /* edge can be on plane or parallel. checking if edge on plane(both points lie on plane)*/
//            if (point_onPlane(v2, nomPL, ptPL) && point_onPlane(v1, nomPL, ptPL))
//            {
//                intersects.push_back(vec1);
//                intersects.push_back(vec2);
//                return intersects;
//            }
//        }
//        else
//        {
//
//            double r = (nomPL.dot(ptPL - v1)) / (nomPL.dot(v2 - v1));
//            if (r >= 0 & r <= 1)
//            {
//
//                Eigen::Vector3d pt = v1 + r * (v2 - v1);
//                intersects.push_back(pt);
//            }
//        }
//    }
//    /* checking if the vertices are duplicate. could be if the triangle is only touching one vertex */
//    if (intersects.size() == 2 && intersects[0] == intersects[1])
//    {
//        intersects.clear();
//    }
//    return intersects;
//}
//
//bool plane_triangleIntersect(Eigen::Vector3d &vec1, Eigen::Vector3d &vec2, Eigen::Vector3d &vec3, Eigen::Vector3d &nomPL, Eigen::Vector3d &ptPL)
//{
//
//    /* check whether the triangle cross the plane or not */
//    double a = (vec1 - ptPL).dot(nomPL);
//    double b = (vec2 - ptPL).dot(nomPL);
//    double c = (vec3 - ptPL).dot(nomPL);
//
//    if (a > 0 && b > 0 && c > 0)
//    { //    /* checking if any two scalars are diff than the other */
//        return false;
//    }
//    else if (a < 0 && b < 0 && c < 0)
//    {
//        return false;
//    }
//    else
//    {
//        return true; // vertices can be on the plane as well
//    }
//}
//
//bool point_onPlane(Eigen::Vector3d &vec1, Eigen::Vector3d &nomPL, Eigen::Vector3d &ptPL)
//{
//    /* checking if point lie on plane. satisfy plane equ */
//    if ((vec1 - ptPL).dot(nomPL) == 0)
//    {
//        return true;
//    }
//    else
//    {
//        return false;
//    }
//}
//
//void updatePoints(std::vector<Eigen::RowVector3d> &outputVertices, Eigen::Vector4d &size)
//{
//    /* rotating the points to make them align to z axis */
//
//    double diff_deg = acos(Eigen::RowVector3d(0, 0, 1).dot(new_normal)); // align to z axis
//    Eigen::RowVector3d axis = new_normal.cross(Eigen::RowVector3d(0, 0, 1));
//
//    Eigen::MatrixXd inputV, outputV;
//    Eigen::Matrix<double, 4, 4> T;
//    Eigen::Matrix<double, 4, 4> R;
//    T.setIdentity();
//
//    inputV.resize(outputVertices.size(), 4);
//    outputV.resize(outputVertices.size(), 4);
//
//    for (int i = 0; i < outputVertices.size(); i++)
//    { // convert to eigen matrix useful to find the centroid
//        inputV.row(i) << outputVertices.at(i)(0), outputVertices.at(i)(1), outputVertices.at(i)(2), 1;
//    }
//
//    Eigen::Vector4d max_contour = inputV.colwise().maxCoeff();
//    Eigen::Vector4d min_contour = inputV.colwise().minCoeff();
//    Eigen::Vector4d centroid = (min_contour - max_contour) / 2;
//
//    T.col(3) << centroid(0), centroid(1), centroid(2), 1;
//
//    MyQuaternion Q(axis, diff_deg);
//    R = Q.Quat_to_Rotmatrix();
//
//    outputV = T * R * T.inverse() * inputV.transpose();
//
//    /* translate so that all vertices are +ve with some buffer*/
//    min_contour = outputV.rowwise().minCoeff();
//    T.col(3) << -min_contour(0) + 50, -min_contour(1) + 50, -min_contour(2) + 50, 1;
//
//    outputV = T * outputV;
//    outputV.transposeInPlace();
//
//    Eigen::RowVector3d newvertices;
//
//    for (int i = 0; i < outputVertices.size(); i++)
//    {
//
//        newvertices << outputV(i, 0), outputV(i, 1), outputV(i, 2);
//        outputVertices.at(i) = newvertices;
//    }
//
//    /* to get the size of the final image */
//    size = outputV.colwise().maxCoeff();
//}
//
//std::vector<Eigen::Vector3d> rayBB_intersection(Eigen::Vector3d &min_BB, Eigen::Vector3d &max_BB, Eigen::Vector3d &ray_direction)
//{
//    /* copmuting the intersection of the new normal with the BB */
//
//    std::vector<Eigen::Vector3d> BB_intersects;
//    Eigen::Vector3d temp;
//
//    if (abs(ray_direction.dot(Eigen::Vector3d(1, 0, 0))) == 1)
//    { // direction parallel to x axis
//
//        BB_intersects.push_back(min_BB);
//
//        temp << max_BB(0), min_BB(1), min_BB(2);
//        BB_intersects.push_back(temp);
//    }
//    else if (abs(ray_direction.dot(Eigen::Vector3d(0, 1, 0))) == 1)
//    { // direction parallel to y axis
//
//        BB_intersects.push_back(min_BB);
//
//        temp << min_BB(0), max_BB(1), min_BB(2);
//        BB_intersects.push_back(temp);
//    }
//    else if (abs(ray_direction.dot(Eigen::Vector3d(0, 0, 1))) == 1)
//    { // direction parallel to z axis
//
//        BB_intersects.push_back(min_BB);
//
//        temp << min_BB(0), min_BB(1), max_BB(2);
//        BB_intersects.push_back(temp);
//    }
//    else
//    {
//
//        double r;
//        double r1 = 1000000.0;
//        Eigen::Vector3d ray_pt1;
//        ray_pt1 = min_BB + ray_direction; // new pt in same direction
//
//        if ((Eigen::Vector3d(1, 0, 0).dot(ray_pt1 - min_BB)) != 0)
//        {
//
//            r = (Eigen::Vector3d(1, 0, 0).dot(max_BB - min_BB)) / (Eigen::Vector3d(1, 0, 0).dot(ray_pt1 - min_BB));
//            if (r1 > r)
//            {
//                r1 = r;
//            }
//        }
//
//        if ((Eigen::Vector3d(0, 1, 0).dot(ray_pt1 - min_BB)) != 0)
//        {
//
//            r = (Eigen::Vector3d(0, 1, 0).dot(max_BB - min_BB)) / (Eigen::Vector3d(0, 1, 0).dot(ray_pt1 - min_BB));
//
//            if (r1 > r)
//            {
//                r1 = r;
//            }
//        }
//
//        if ((Eigen::Vector3d(0, 0, 1).dot(ray_pt1 - min_BB)) != 0)
//        {
//
//            r = (Eigen::Vector3d(0, 0, 1).dot(max_BB - min_BB)) / (Eigen::Vector3d(0, 0, 1).dot(ray_pt1 - min_BB));
//
//            if (r1 > r)
//            {
//                r1 = r;
//            }
//        }
//        BB_intersects.push_back(min_BB);
//        BB_intersects.push_back(min_BB + r1 * ray_direction);
//    }
//    return BB_intersects;
//}
//
//void writePLY(std::vector<Eigen::RowVector3d> VERTICES, std::vector<Eigen::RowVector3d> COLOR, std::vector<Eigen::RowVector2d> EDGES, std::string directory)
//{
//
//    FILE *PLYfile = fopen(directory.c_str(), "w");
//    if (!PLYfile)
//    {
//        std::cout << "open file failed!" << std::endl;
//    }
//
//    fprintf(PLYfile, "ply\n");
//    fprintf(PLYfile, "format ascii 1.0\n");
//    fprintf(PLYfile, "element vertex %d\n", (int)VERTICES.size());
//    fprintf(PLYfile, "property float x\n");
//    fprintf(PLYfile, "property float y\n");
//    fprintf(PLYfile, "property float z\n");
//    fprintf(PLYfile, "property uchar red\n");
//    fprintf(PLYfile, "property uchar green\n");
//    fprintf(PLYfile, "property uchar blue\n");
//    fprintf(PLYfile, "element faces %d\n", (int)EDGES.size());
//    fprintf(PLYfile, "property list uchar int vertex_indices\n");
//    fprintf(PLYfile, "end_header\n");
//
//    for (int i = 0; i < VERTICES.size(); i++)
//    {
//
//        fprintf(PLYfile, "%.3f %.3f %.3f %3.0f %3.0f %3.0f\n", VERTICES.at(i)(0), VERTICES.at(i)(1), VERTICES.at(i)(2), COLOR.at(i)(0) * 255, COLOR.at(i)(1) * 255, COLOR.at(i)(2) * 255);
//    }
//    for (int i = 0; i < EDGES.size(); i++)
//    {
//
//        fprintf(PLYfile, "2 %d %d\n", (int)EDGES.at(i)(0), (int)EDGES.at(i)(1));
//    }
//    fclose(PLYfile);
//}
//
//void writePNG_FILE(std::vector<Eigen::RowVector3d> VERTICES, std::vector<Eigen::RowVector3d> COLORS, std::vector<Eigen::RowVector2d> EDGES, std::string PNG_FILE)
//{
//
//    // png file configure
//    const char *png_file = PNG_FILE.c_str();
//    int width = 640;
//    int height = 480;
//
//    std::cout << "VERTICES size : " << VERTICES.size() << std::endl;
//    std::cout << "COLORS size : " << COLORS.size() << std::endl;
//    std::cout << "EDGES size : " << EDGES.size() << std::endl;
//
//    Eigen::RowVector3d last_color = COLORS[0];
//    std::vector<int> partitions;
//
//    for(int i = 0; i < COLORS.size(); ++i) {
//        if(COLORS[i] != last_color) {
//            partitions.emplace_back(i);
//        }
//        last_color = COLORS[i];
//    }
//    std::cout << "partitions : " << "0~" << partitions[0] - 1 << std::endl;
//
//    Eigen::MatrixXd V1(partitions[0], 2);
//    Eigen::MatrixXi E1(partitions[0], 2);
//    V1.setZero();
//    E1.setZero();
//
//    for (int i = 0; i < partitions[0]; ++i)
//    {
//        V1(i,0) = VERTICES.at(i)(0);
//        V1(i,1) = VERTICES.at(i)(2);
//    }
//    for (int i = 0; i < partitions[0]; ++i)
//    {
//        E1(i,0) = i;
//        E1(i,1) = i+1;
//    }
//    E1(partitions[0] - 1, 1) = E1(0,0);
//    std::cout << "V1:\n" << V1 << std::endl;
//    std::cout << "E1:\n" << E1 << std::endl;
//
//    Eigen::MatrixXd V2;
//    Eigen::MatrixXi F;
//
//    igl::triangle::triangulate(V1, E1, Eigen::MatrixXd(), "a1.000q", V2, F);
//    Eigen::MatrixXd V(V2.rows(), 3);
//    V.setZero();
//    V.col(0) = V2.col(0);
//    V.col(1) = V2.col(1);
//    std::cout << V << std::endl;
//
//     Eigen::VectorXd K;
//     // Compute integral of Gaussian curvature
//     igl::gaussian_curvature(V, F, K);
//     // Compute mass matrix
//     Eigen::SparseMatrix<double> M, Minv;
//     igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
//     igl::invert_diag(M, Minv);
//     // Divide by area to get integral average
//     K = (Minv * K).eval();
//
//     // embree object
//     igl::embree::EmbreeRenderer er;
//     er.set_mesh(V, F, true);
//
//     er.set_data(K, igl::COLOR_MAP_TYPE_JET);
//     er.set_colors(Eigen::RowVector3d(0.0, 1.0, 1.0));
//
//     Eigen::Matrix3d rot_matrix;
//
//     // specify rotation
//     rot_matrix = Eigen::AngleAxisd(0 * igl::PI / 180.0, Eigen::Vector3d::UnitX()) * Eigen::AngleAxisd(0 * igl::PI / 180.0, Eigen::Vector3d::UnitY()) * Eigen::AngleAxisd(0 * igl::PI / 180.0, Eigen::Vector3d::UnitZ());
//     er.set_rot(rot_matrix);
//
//     er.set_zoom(1.5);
//     er.set_orthographic(true);
//
//     Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> R(width, height);
//     Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> G(width, height);
//     Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> B(width, height);
//     Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> A(width, height);
//
//     // render view using embree
//     er.render_buffer(R, G, B, A);
//
//     std::cout << "Rendered scene saved to " << png_file << std::endl;
//
//     // save to PNG file
//     igl::png::writePNG(R, G, B, A, png_file);
//}


#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <igl/readOFF.h>

#include <igl/lexicographic_triangulation.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/triangle/triangulate.h>
#include <igl/gaussian_curvature.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/read_triangle_mesh.h>
#include <igl/png/writePNG.h>
#include <igl/PI.h>
#include <Eigen/Geometry>
// embree
#include <igl/embree/EmbreeRenderer.h>

#include <igl/predicates/predicates.h>
#include <igl/delaunay_triangulation.h>

#include <iostream>

int orient2dPredicates(const double* pa, const double* pb, const double* pc)
{
    const Eigen::Vector2d a(pa[0], pa[1]);
    const Eigen::Vector2d b(pb[0], pb[1]);
    const Eigen::Vector2d c(pc[0], pc[1]);

    const auto result = igl::predicates::orient2d<Eigen::Vector2d>(a, b, c);

    if (result == igl::predicates::Orientation::POSITIVE) {
        return 1;
    }
    else if (result == igl::predicates::Orientation::NEGATIVE) {
        return -1;
    }
    else {
        return 0;
    }
}

int inCirclePredicates(const double* pa, const double* pb, const double* pc, const double* pd)
{
    const Eigen::Vector2d a(pa[0], pa[1]);
    const Eigen::Vector2d b(pb[0], pb[1]);
    const Eigen::Vector2d c(pc[0], pc[1]);
    const Eigen::Vector2d d(pd[0], pd[1]);

    const auto result = igl::predicates::incircle(a, b, c, d);

    if (result == igl::predicates::Orientation::INSIDE) {
        return 1;
    }
    else if (result == igl::predicates::Orientation::OUTSIDE) {
        return -1;
    }
    else {
        return 0;
    }
}

int main() {


    // How I call the delaunay triangulation function
    Eigen::MatrixXd vertices(6, 2);
    Eigen::MatrixXi edges(6, 2);

    
    // both centered at origin
    vertices << -2,0, -2,2, 2,2, 2,-2, 0,-2, 0,0;
    edges << 0,1, 1,2, 2,3, 3,4, 4,5, 5,0;
        
    Eigen::MatrixXi faces;
    Eigen::MatrixXi vertices2;
    igl::delaunay_triangulation(vertices, orient2dPredicates, inCirclePredicates, faces);
    //igl::triangle::triangulate(vertices, edges, Eigen::MatrixXd(), "a0.005q", vertices2, faces);

    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(vertices, faces);
    viewer.launch();
}