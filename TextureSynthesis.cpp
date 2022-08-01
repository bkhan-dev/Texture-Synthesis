#define MAXIMUM_LINE_LENGTH 1024
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

#include "Cartesian3.h"

// will store information for a pixel in the image
struct imagePixel {
    Cartesian3 coords; // coordinates in the image
    Cartesian3 col; // colour of vertex
    Cartesian3 n; // normal
    imagePixel(){coords = Cartesian3(); col = Cartesian3(); n = Cartesian3();}
    imagePixel(Cartesian3 p_coords, Cartesian3 p_col, Cartesian3 p_n)
    {
        coords = p_coords; // convert to image coordinates on creation
        col = p_col;
        n = p_n;
    }
};

// reads in the input file - using code from the renderer
bool readModel(std::istream &modelStream);
// will rasterise a triangle in the image
void rasterise(imagePixel v1, imagePixel v2, imagePixel v3);

// will hold the texture from the model
Cartesian3 image[1024][1024];

// will hold the normap map
Cartesian3 normalMap[1024][1024];

// storing model data
std::vector<Cartesian3> vertices;
std::vector<Cartesian3> colours;
std::vector<Cartesian3> normals;
std::vector<Cartesian3> textureCoords;

std::vector<unsigned int> faceVertices;
std::vector<unsigned int> faceColours;
std::vector<unsigned int> faceNormals;
std::vector<unsigned int> faceTexCoords;

int main(int argc, char *argv[])
{
    // Exit program if no model is given, or if more than 1 model is given
    if(argc != 2){
        std::cout << "Usage: " << argv[0] << "  PATH_TO/<modelname>.obj\n";
        return 0;
    };

    // Open given .obj file
    std::ifstream model(argv[1]);

    // Exit if the file could not be opened
    if(!model.is_open()){
        std::cout << "Could not open file.\n";
        return -1;
    }

    // read in model data using code from renderer
    readModel(model);

    std::cout << "\nGenerating texture and normal map..." << std::endl;

    // go over every face
    for(int i = 0; i < faceVertices.size(); i+=3){

        // get the pixel info - conversions to proper number ranges also happen here
        imagePixel v1(textureCoords[faceTexCoords[i]] * 1023, 255.0f * colours[faceColours[i]], 255.0f * ((normals[faceNormals[i]] + Cartesian3(1.0f,1.0f,1.0f)) / 2.0f));
        imagePixel v2(textureCoords[faceTexCoords[i+1]] * 1023, 255.0f * colours[faceColours[i+1]], 255.0f * ((normals[faceNormals[i+1]] + Cartesian3(1.0f,1.0f,1.0f)) / 2.0f));
        imagePixel v3(textureCoords[faceTexCoords[i+2]] * 1023, 255.0f * colours[faceColours[i+2]], 255.0f * ((normals[faceNormals[i+2]] + Cartesian3(1.0f,1.0f,1.0f)) / 2.0f));

        // rasterise the triangle in the images
        rasterise(v1,v2,v3);
    }

    std::cout << "Finished generating." << std::endl;


    std::cout << "\nCreating PPM files" << std::endl;


    // write out ppm for the texture
    std::ofstream ppm("texture.ppm");

    if(!ppm.is_open()){
        return -1;
    }

    ppm << "P3" << std::endl;
    ppm << 1024 << " " << 1024 << std::endl;
    ppm << 255 << std::endl;

    for(int row = 1023; row >= 0; row--)
        for(int col = 0; col < 1024; col++)
            ppm << image[row][col].x << " " << image[row][col].y << " " << image[row][col].z << "  ";

    ppm.close();


    // write out ppm for the normal map
    std::ofstream ppmN("normalMap.ppm");

    if(!ppmN.is_open()){
        return -1;
    }

    ppmN << "P3" << std::endl;
    ppmN << 1024 << " " << 1024 << std::endl;
    ppmN << 255 << std::endl;

    for(int row = 1023; row >= 0; row--)
        for(int col = 0; col < 1024; col++)
            ppmN << normalMap[row][col].x << " " << normalMap[row][col].y << " " << normalMap[row][col].z << "  ";

    ppmN.close();

    std::cout << "Images created.\n" << std::endl;

    return 0;
}


void rasterise(imagePixel v0i, imagePixel v1i, imagePixel v2i){
    // get points
    Cartesian3 v0((int)v0i.coords.x, (int)v0i.coords.y, 0);
    Cartesian3 v1((int)v1i.coords.x, (int)v1i.coords.y, 0);
    Cartesian3 v2((int)v2i.coords.x, (int)v2i.coords.y, 0);

    // set the bounding box to coords of v0 as a starting point
    int minX, maxX, minY, maxY;
    minX = maxX = v0.x;
    minY = maxY = v0.y;

    // test against other vertices
    // minimum X
    if(v1.x < minX) minX = v1.x;
    if(v2.x < minX) minX = v2.x;
    // minimum y
    if(v1.y < minY) minY = v1.y;
    if(v2.y < minY) minY = v2.y;
    // maximum X
    if(v1.x > maxX) maxX = v1.x;
    if(v2.x > maxX) maxX = v2.x;
    // maximum Y
    if(v1.y > maxY) maxY = v1.y;
    if(v2.y > maxY) maxY = v2.y;

    // compute information about the triangle
    Cartesian3 v0v1 = v1 - v0;
    Cartesian3 v0v2 = v2 - v0;
    Cartesian3 normal = v0v1.cross(v0v2);
    float areaABC = normal.dot(v0v1.cross(v0v2));

    // go over every pixel in the bounding box
    for(int x = minX; x < maxX; x++){
        for(int y = minY; y < maxY; y++){
            // get current point
            Cartesian3 point(x,y,0);

            // Get area of two sub triangles
            float areaPBC = normal.dot((v1 - point).cross(v2 - point));
            float areaPCA = normal.dot((v2 - point).cross(v0 - point));

            // get alpha beta gamma
            float alpha = areaPBC / areaABC;
            float beta = areaPCA / areaABC;
            float gamma = 1.0f - alpha - beta;

            // check if it is inside the triangle or not
            if ((alpha < 0.0) || (beta < 0.0) || (gamma < 0.0))
                continue;

            // get the weighted colour, this is for the texture image
            Cartesian3 col = (alpha * v0i.col) + (beta * v1i.col) + (gamma * v2i.col);
            // convert to int
            col.x = (int)col.x;
            col.y = (int)col.y;
            col.z = (int)col.z;

            // set the colour in the texture image
            image[y][x] = col;

            // get the weighted colour of the normal which has already been mapped to RGB
            Cartesian3 nCol((alpha * v0i.n) + (beta * v1i.n) + (gamma * v2i.n));
            // convert to int
            nCol.x = (int)nCol.x;
            nCol.y = (int)nCol.y;
            nCol.z = (int)nCol.z;

            // set the colour in the normal map
            normalMap[y][x] = nCol;

        }
    }

}

bool readModel(std::istream &modelStream)
{
     // create a read buffer
    char readBuffer[MAXIMUM_LINE_LENGTH];

    // the rest of this is a loop reading lines & adding them in appropriate places
    while (true)
        { // not eof
        // character to read
        char firstChar = modelStream.get();

//         std::cout << "Read: " << firstChar << std::endl;

        // check for eof() in case we've run out
        if (modelStream.eof())
            break;

        // otherwise, switch on the character we read
        switch (firstChar)
            { // switch on first character
            case '#':       // comment line
                // read and discard the line
                modelStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
                break;

            case 'v':       // vertex data of some type
                { // some sort of vertex data
                // retrieve another character
                char secondChar = modelStream.get();

                // bail if we ran out of file
                if (modelStream.eof())
                    break;

                // now use the second character to choose branch
                switch (secondChar)
                    { // switch on second character
                    case ' ':       // space - indicates a vertex
                        { // vertex read
                        Cartesian3 vertex;
                        modelStream >> vertex;
                        vertices.push_back(vertex);
//                         std::cout << "Vertex " << vertex << std::endl;
                        break;
                        } // vertex read
                    case 'c':       // c indicates colour
                        { // normal read
                        Cartesian3 colour;
                        modelStream >> colour;
                        colours.push_back(colour);
//                         std::cout << "Colour " << colour << std::endl;
                        break;
                        } // normal read
                    case 'n':       // n indicates normal vector
                        { // normal read
                        Cartesian3 normal;
                        modelStream >> normal;
                        normals.push_back(normal);
//                         std::cout << "Normal " << normal << std::endl;
                        break;
                        } // normal read
                    case 't':       // t indicates texture coords
                        { // tex coord
                        Cartesian3 texCoord;
                        modelStream >> texCoord;
                        textureCoords.push_back(texCoord);
//                         std::cout << "Tex Coords " << texCoord << std::endl;
                        break;
                        } // tex coord
                    default:
                        break;
                    } // switch on second character
                break;
                } // some sort of vertex data

            case 'f':       // face data
                { // face
				// make a hard assumption that we have a single triangle per line
                unsigned int vertexID;
				unsigned int colourID;
                unsigned int normalID;
				unsigned int texCoordID;

                // read in three vertices
				for (unsigned int vertex = 0; vertex < 3; vertex++)
					{ // per vertex
					// read a vertex ID
					modelStream >> vertexID;
					// read and discard the slash
					modelStream.get();
					// read a colour ID
					modelStream >> colourID;
					// read and discard the slash
					modelStream.get();
					// read a vertex ID
					modelStream >> texCoordID;
					// read and discard the slash
					modelStream.get();
					// read a vertex ID
					modelStream >> normalID;

// 					std::cout << "Face " << vertexID << "/" << colourID << "/" << texCoordID << "/" << normalID << std::endl;

					// subtract one and store them (OBJ uses 1-based numbering)
					faceVertices.push_back(vertexID-1);
					faceColours.push_back(colourID-1);
					faceNormals.push_back(normalID-1);
					faceTexCoords.push_back(texCoordID-1);
					} // per vertex
				break;
                } // face

            // default processing: do nothing
            default:
                break;

            } // switch on first character

        } // not eof
}
