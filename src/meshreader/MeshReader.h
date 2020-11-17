/**
 * @file
 *  This file is part of PUMGen
 *
 *  For conditions of distribution and use, please see the copyright
 *  notice in the file 'COPYING' at the root directory of this package
 *  and the copyright notice at https://github.com/SeisSol/PUMGen
 *
 * @copyright 2017 Technical University of Munich
 * @author Sebastian Rettenberger <sebastian.rettenberger@tum.de>
 */

#ifndef MESH_READER_H
#define MESH_READER_H

#include <fstream>

#include "utils/logger.h"

class MeshReader {
protected:
    struct FileSection {
        /** Start of the section */
        size_t seekPosition;
        /** Number of elements (lines) in the section */
        unsigned int nLines;
        /** Line size */
        size_t lineSize;
    };

    std::ifstream m_mesh;

public:
    MeshReader() {}

    MeshReader(const char* meshFile) { open(meshFile); }

    virtual ~MeshReader() {}

    virtual void open(const char* meshFile) {
        m_mesh.open(meshFile);

        // Files are readable?
        if (!m_mesh)
            logError() << "Could not open mesh file" << meshFile;
    }

    virtual unsigned int nVertices() const = 0;
    virtual unsigned int nElements() const = 0;
    /**
     * @return Number of boundary faces
     */
    virtual unsigned int nBoundaries() const = 0;

    /**
     * Reads vertices from start tp start+count from the file and stores them in
     * vertices
     *
     * @param vertices Buffer to store coordinates of each vertex. The caller is
     * responsible for allocating the buffer. The size of the buffer must be
     * count*dimensions
     */
    virtual void readVertices(unsigned int start, unsigned int count, double* vertices) = 0;

    /**
     * Read all vertices
     *
     * @see readVertices(unsigned int, unsigned int, double*)
     */
    void readVertices(double* vertices) {
        logInfo() << "Reading vertices";
        readVertices(0, nVertices(), vertices);
    }

    /**
     * Reads elements from start to start+count from the file and stores them in
     * elements
     *
     * @param elements Buffer to store vertices of each element. The caller is
     * responsible for allocating the buffer. The Size of the buffer must be
     * count*vertices_per_element.
     */
    virtual void readElements(unsigned int start, unsigned int count, int* elements) = 0;

    /**
     * Reads all elements
     *
     * @see readElements(unsigned int, unsigned int, unsinged int*)
     */
    void readElements(int* elements) {
        logInfo() << "Reading elements";
        readElements(0, nElements(), elements);
    }
};

#endif // MESH_READER_H
