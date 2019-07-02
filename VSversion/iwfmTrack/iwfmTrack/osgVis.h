#pragma once

#include <osg/Geometry>
#include <osg/Geode>
#include <osgUtil/SmoothingVisitor>

#include "iwfmTree.h"


class osgMesh {
public:
	osgMesh(iwfmMesh& iwfmMSH);
	osg::Geode* getScene();

private:
	osg::ref_ptr<osg::Geode> root = new osg::Geode;
};


osgMesh::osgMesh(iwfmMesh& M) {
	osg::ref_ptr<osg::Vec3Array> vert_top = new osg::Vec3Array;
	osg::ref_ptr<osg::Vec3Array> vert_bot = new osg::Vec3Array;

	int l = M.Nlay() - 1;

	for (int i = 0; i < M.Nnodes(); ++i) {
		vert_top->push_back(osg::Vec3(static_cast<float>(M.getXY(i, 0)),
									  static_cast<float>(M.getXY(i, 1)),
									  static_cast<float>(M.getZ(i, 0))));
		vert_bot->push_back(osg::Vec3(static_cast<float>(M.getXY(i, 0)),
									  static_cast<float>(M.getXY(i, 1)),
									  static_cast<float>(M.getZ(i, l))));
	}

	osg::ref_ptr<osg::DrawElementsUInt> indices = new osg::DrawElementsUInt(GL_TRIANGLES);
	for (int i = 0; i < M.Nelem(); ++i) {
		if (M.getMSH(i, 3) == 0) {
			(*indices).push_back(M.getMSH(i, 0) - 1);
			(*indices).push_back(M.getMSH(i, 1) - 1);
			(*indices).push_back(M.getMSH(i, 2) - 1);
		}
		else {
			(*indices).push_back(M.getMSH(i, 0) - 1);
			(*indices).push_back(M.getMSH(i, 1) - 1);
			(*indices).push_back(M.getMSH(i, 2) - 1);
			(*indices).push_back(M.getMSH(i, 2) - 1);
			(*indices).push_back(M.getMSH(i, 3) - 1);
			(*indices).push_back(M.getMSH(i, 0) - 1);
		}
	}

	osg::ref_ptr<osg::Geometry> top_mesh = new osg::Geometry;
	osg::ref_ptr<osg::Geometry> bot_mesh = new osg::Geometry;
	top_mesh->setVertexArray(vert_top.get());
	top_mesh->addPrimitiveSet(indices.get());
	bot_mesh->setVertexArray(vert_bot.get());
	bot_mesh->addPrimitiveSet(indices.get());

	osgUtil::SmoothingVisitor::smooth(*top_mesh);

	root->addDrawable(top_mesh.get());
	root->addDrawable(bot_mesh.get());
}

osg::Geode* osgMesh::getScene() {
	return root.release();
}