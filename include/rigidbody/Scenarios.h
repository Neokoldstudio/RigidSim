#pragma once

#include "rigidbody/RigidBody.h"
#include "rigidbody/RigidBodySystem.h"
#include "joint/Spherical.h"
#include "joint/Distance.h"
#include "util/Types.h"
#include "util/MeshAssets.h"

#include <cstdlib>
#include <Eigen/Dense>


class Scenarios
{
public:
    // Box filled with balls.
    //
    static void createMarbleBox(RigidBodySystem& rigidBodySystem)
    {
        rigidBodySystem.clear();
        polyscope::removeAllStructures();

        std::cout << "Loading marble box scenario" << std::endl;

        // Create two layers of "marbles", in a grid layout.
        //
        for (int i = 0; i < 9; ++i)
        {
            for (int j = 0; j < 9; ++j)
            {
                RigidBody* body1 = new RigidBody(1.0f, new Sphere(0.5f), createSphere(0.5f));
                body1->x = { -4.0f + (float)i * 1.0f, 2.0f, -4.0f + (float)j * 1.0f };
                body1->xdot = Eigen::Vector3f::Random();
                rigidBodySystem.addBody(body1);
                body1->mesh->setSurfaceColor({ 1.0f, 0.1f, 0.1f });
                body1->mesh->setTransparency(0.6f);
                RigidBody* body2 = new RigidBody(1.0f, new Sphere(0.5f), createSphere(0.5f));
                body2->x = { -4.0f + (float)i * 1.0f, 3.0f, -4.0f + (float)j * 1.0f };
                body2->xdot = Eigen::Vector3f::Random();
                rigidBodySystem.addBody(body2);
                body2->mesh->setSurfaceColor({ 1.0f, 0.1f, 0.1f });
                body2->mesh->setTransparency(0.6f);
            }
        }

        // Create the box to hold the marbles.
        RigidBody* body0 = new RigidBody(1.0f, new Box(Eigen::Vector3f(0.4f, 4.0f, 10.0f)), createBox(Eigen::Vector3f(0.4f, 4.0f, 10.0f)));
        RigidBody* body1 = new RigidBody(1.0f, new Box(Eigen::Vector3f(0.4f, 4.0f, 10.0f)), createBox(Eigen::Vector3f(0.4f, 4.0f, 10.0f)));
        RigidBody* body2 = new RigidBody(1.0f, new Box(Eigen::Vector3f(0.4f, 4.0f, 10.0f)), createBox(Eigen::Vector3f(0.4f, 4.0f, 10.0f)));
        RigidBody* body3 = new RigidBody(1.0f, new Box(Eigen::Vector3f(0.4f, 4.0f, 10.4f)), createBox(Eigen::Vector3f(0.4f, 4.0f, 10.0f)));
        RigidBody* body4 = new RigidBody(1.0f, new Box(Eigen::Vector3f(10.0f, 0.4f, 10.0f)), createBox(Eigen::Vector3f(10.0f, 0.4f, 10.0f)));
        body0->fixed = true;
        body1->fixed = true;
        body2->fixed = true;
        body3->fixed = true;
        body4->fixed = true;
        body0->mesh->setSurfaceColor({ 0.6f, 0.6f, 0.6f })->setSmoothShade(false)->setTransparency(0.4f);
        body1->mesh->setSurfaceColor({ 0.6f, 0.6f, 0.6f })->setSmoothShade(false)->setTransparency(0.4f);
        body2->mesh->setSurfaceColor({ 0.6f, 0.6f, 0.6f })->setSmoothShade(false)->setTransparency(0.4f);
        body3->mesh->setSurfaceColor({ 0.6f, 0.6f, 0.6f })->setSmoothShade(false)->setTransparency(0.4f);
        body4->mesh->setSurfaceColor({ 0.6f, 0.6f, 0.6f })->setSmoothShade(false)->setTransparency(0.4f);
        body0->x = { 4.75f, 2.0f, 0.0f };
        body1->x = { -4.75f, 2.0f, 0.0f };
        body2->x = { 0.0f, 2.0f, 4.75f };
        body2->q = Eigen::AngleAxisf(1.57f, Eigen::Vector3f(0, 1, 0));
        body3->x = { 0.0f, 2.0f, -4.75f };
        body3->q = Eigen::AngleAxisf(1.57f, Eigen::Vector3f(0, 1, 0));
        body4->x = { 0.0f, 0.0f, 0.0f };

        rigidBodySystem.addBody(body0);
        rigidBodySystem.addBody(body1);
        rigidBodySystem.addBody(body2);
        rigidBodySystem.addBody(body3);
        rigidBodySystem.addBody(body4);
    }

    // Sphere-sphere distance test
    //
    static void createSphereSphereDistance(RigidBodySystem& rigidBodySystem)
    {
        rigidBodySystem.clear();
        polyscope::removeAllStructures();

        std::cout << "Loading sphere-sphere distance test." << std::endl;

        // Create a sphere.
        RigidBody* sphere1 = new RigidBody(1.0f, new Sphere(0.5f), createSphere(0.5f));
        sphere1->x = { 0.0f, 4.0f, 0.0f };
        sphere1->fixed = true;

        // Create another sphere.
        RigidBody* sphere2 = new RigidBody(1.0f, new Sphere(0.5f), createSphere(0.5f));
        sphere2->x = { -2.0f, 4.0f, 0.0f };

        rigidBodySystem.addBody(sphere1);
        rigidBodySystem.addBody(sphere2);

        Distance* joint = new Distance(sphere1, sphere2, { 0.0f, -0.5f, 0.0f }, { 0.0f, -0.5f, 0.0f }, 2.0f);
        rigidBodySystem.addJoint(joint);
    }

    // Sphere-sphere distance test
    //
    static void createSphereInsideBox(RigidBodySystem& rigidBodySystem)
    {
        rigidBodySystem.clear();
        polyscope::removeAllStructures();

        std::cout << "Loading sphere-inside-box test." << std::endl;

        // Create a box.
        const Eigen::Vector3f dim(2.0f, 2.0f, 2.0f);
        RigidBody* box = new RigidBody(1.0f, new Box(dim), createBox(dim));
        box->x = { 0.0f, 2.0f, 0.0f };
        box->fixed = true;
        box->mesh->setTransparency(0.4f);

        // Create another sphere.
        RigidBody* sphere = new RigidBody(1.0f, new Sphere(0.5f), createSphere(0.5f));
        sphere->x = { 0.0f, 2.5f, 0.0f };
        sphere->mesh->setTransparency(0.6f);

        rigidBodySystem.addBody(box);
        rigidBodySystem.addBody(sphere);
    }


    // Simple sphere falling on a box.
    //
    static void createSphereOnBox(RigidBodySystem& rigidBodySystem)
    {
        rigidBodySystem.clear();
        polyscope::removeAllStructures();

        std::cout << "Loading sphere-on-box scenario." << std::endl;

        // Create a sphere.
        RigidBody* bodySphere = new RigidBody(1.0f, new Sphere(0.5f), createSphere(0.5f));
        bodySphere->x.y() = 4.0f;
        bodySphere->omega = Eigen::Vector3f(10.0f, 0.0f, 0.0f);
        bodySphere->mesh->setTransparency(0.8f);

        // Create a box that will act as the ground.
        RigidBody* bodyBox = new RigidBody(1.0f, new Box(Eigen::Vector3f(10.0f, 0.4f, 10.0f)), createBox(Eigen::Vector3f(10.0f, 0.4f, 10.0f)));
        bodyBox->fixed = true;

        rigidBodySystem.addBody(bodySphere);
        rigidBodySystem.addBody(bodyBox);


        bodySphere->mesh->setSurfaceColor({ 0.1f, 1.0f, 0.2f })->setEdgeWidth(1.0f)->setTransparency(0.6f);
        bodyBox->mesh->setSurfaceColor({ 0.2f, 0.2f, 0.2f })->setSmoothShade(false)->setTransparency(0.4f);
    }

    // Box falling on a plane.
    //
    static void createBoxOnPlane(RigidBodySystem& rigidBodySystem)
    {
        rigidBodySystem.clear();
        polyscope::removeAllStructures();

        std::cout << "Loading box-on-plane scenario." << std::endl;

        // Create a box.
        const Eigen::Vector3f dim(1.0f, 1.0f, 1.0f);
        RigidBody* bodyBox = new RigidBody(1.0f, new Box(dim), createBox(dim));
        bodyBox->x = { -1.0f, 2.0f, 0.0f };
        bodyBox->q = Eigen::AngleAxisf(-0.5236, Eigen::Vector3f(0, 0, 1));

        // Create a plane that will act as the ground.
        const Eigen::Vector3f n({ 0.5f, 0.866f, 0.0f });
        const Eigen::Vector3f p({ 0.0f, 0.0f, 0.0f });
        RigidBody* plane = new RigidBody(1.0f, new Plane(p, n), createPlane(p, n));
        plane->fixed = true;

        rigidBodySystem.addBody(bodyBox);
        rigidBodySystem.addBody(plane);

        bodyBox->mesh->setSurfaceColor({ 0.1f, 1.0f, 0.2f })->setEdgeWidth(1.0f)->setTransparency(0.6f);
        plane->mesh->setSurfaceColor({ 0.2f, 0.2f, 0.2f })->setSmoothShade(false)->setTransparency(0.4f);
    }

    // Cylinder falling on sphere.
    //
    static void createCylinderSphereTest(RigidBodySystem& rigidBodySystem, const Eigen::AngleAxisf& aa)
    {
        rigidBodySystem.clear();
        polyscope::removeAllStructures();

        std::cout << "Loading cylinder on sphere scenario." << std::endl;

        // Create a sphere.
        RigidBody* bodySphere = new RigidBody(1.0f, new Sphere(0.5f), createSphere(0.5f));
        bodySphere->x = { 0.0f, 4.0f, 0.0f };
        bodySphere->omega = { 0.0f, 0.0f, 1.0f };
        bodySphere->mesh->setTransparency(0.8f);

        // Create a vertical cylinder.
        RigidBody* bodyCyl = new RigidBody(1.0f, new Cylinder(2.0f, 0.5f), createCylinder(16, 0.5f, 2.0f) );
        bodyCyl->x = { 0.0f, 1.0f, 0.0f };
        bodyCyl->q = Eigen::Quaternionf(aa);
        bodyCyl->fixed = true;

        rigidBodySystem.addBody(bodySphere);
        rigidBodySystem.addBody(bodyCyl);

        bodySphere->mesh->setSurfaceColor({ 0.1f, 1.0f, 0.2f })->setEdgeWidth(1.0f)->setTransparency(0.4f);
        bodyCyl->mesh->setSurfaceColor({ 0.1f, 0.1f, 1.0f })->setSmoothShade(false)->setTransparency(0.4f);
    }

    // Chain of swinging boxes
    //
    static void createRopeLadder(RigidBodySystem& rigidBodySystem)
    {
        rigidBodySystem.clear();
        polyscope::removeAllStructures();

        std::cout << "Loading rope ladder scenario." << std::endl;

        const int N = 8;

        const Eigen::Vector3f dim = { 1.0f, 0.5f, 0.6f };
        const float dy = 0.8f;

        // Create a box.
        RigidBody* topBox = new RigidBody(1.0f, new Box(dim), createBox(dim));
        topBox->x = { 0.0f, dy * (float)N, 0.0f };
        topBox->fixed = true;
        rigidBodySystem.addBody(topBox);

        RigidBody* parent = topBox;
        for (int i = 0; i < N - 1; ++i)
        {
            // Create the next box in the chain.
            RigidBody* nextBox = nullptr;
            nextBox = new RigidBody(1.0f, new Box(dim), createBox(dim));
            nextBox->x = parent->x - Eigen::Vector3f(0.0f, dy, 0.0f);
            rigidBodySystem.addBody(nextBox);

            // Add a distance constraint between parent->nextBox
            rigidBodySystem.addJoint(new Distance(parent, nextBox, { 0.5f, 0.0f, 0.0f }, { 0.5f, 0.0f, 0.0f }, dy));
            rigidBodySystem.addJoint(new Distance(parent, nextBox, { -0.5f, 0.0f, 0.0f }, { -0.5f, 0.0f, 0.0f }, dy));

            // Update next body
            parent = nextBox;
        }

        parent->xdot = { 0.0f, 0.0, 5.0f };
    }

    // TODO Custom Scenario.
    //
    static void createCustomScenario(RigidBodySystem& rigidBodySystem)
    {
        rigidBodySystem.clear();
        polyscope::removeAllStructures();

        std::cout << "Loading custom scenario." << std::endl;

        // TODO Objective #7 
        // Create your own interesting scenario. Be creative!
        // Your scenario should be more than just a trivial extension 
        // of an existing example.
        //

    }

};
