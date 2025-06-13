#include "viewer/SimViewer.h"

#include "polyscope/curve_network.h"
#include "polyscope/pick.h"
#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/view.h"
#include "imgui.h"

#include <chrono>
#include <iostream>
#include <functional>

#include "contact/Contact.h"
#include "rigidbody/RigidBodySystem.h"
#include "rigidbody/RigidBodyState.h"
#include "rigidbody/Scenarios.h"

using namespace std;

namespace
{
    static RigidBodySystem* m_rigidBodySystem = new RigidBodySystem;

    static void updateRigidBodyMeshes(RigidBodySystem& _rigidBodySystem)
    {
        auto& bodies = _rigidBodySystem.getBodies();
        for(unsigned int k = 0; k < bodies.size(); ++k)
        { 
            if (!bodies[k]->mesh) continue;

            Eigen::Isometry3f tm = Eigen::Isometry3f::Identity();
        
            // copy rotation part
            tm.linear() = bodies[k]->q.toRotationMatrix();

            // copy translation part
            tm.translation() = bodies[k]->x;

            bodies[k]->mesh->setTransform(glm::make_mat4x4(tm.data()));
        }
    }

    static void updateContactPoints(RigidBodySystem& _rigidBodySystem)
    {
        const auto& contacts = _rigidBodySystem.getContacts();
        const unsigned int numContacts = contacts.size();

        if (numContacts == 0)
        {
            polyscope::removePointCloud("contacts");
        }
        else
        {
            Eigen::MatrixXf contactP(numContacts, 3);
            Eigen::MatrixXf contactN(numContacts, 3);

            for (unsigned int i = 0; i < numContacts; ++i)
            {
                contactP.row(i)(0) = contacts[i].p(0); contactP.row(i)(1) = contacts[i].p(1); contactP.row(i)(2) = contacts[i].p(2);
                contactN.row(i)(0) = contacts[i].n(0); contactN.row(i)(1) = contacts[i].n(1); contactN.row(i)(2) = contacts[i].n(2);
            }

            auto pointCloud = polyscope::registerPointCloud("contacts", contactP);

            pointCloud->setPointColor({ 1.0f, 0.0f, 0.0f });
            pointCloud->setPointRadius(0.005);
            pointCloud->addVectorQuantity("normal", contactN)->setVectorColor({ 1.0f, 1.0f, 0.0f })->setVectorLengthScale(0.05f)->setEnabled(true);
        }
    }

    static void updateJointViz(RigidBodySystem& _rigidBodySystem)
    {
        const auto& joints = _rigidBodySystem.getJoints();
        const unsigned int numJoints = joints.size();

        if (numJoints == 0)
        {
            polyscope::removePointCloud("jointsPoint");
            polyscope::removeCurveNetwork("jointsCurve");
        }
        else
        {
            Eigen::MatrixXf jointP(2 * numJoints, 3);
            Eigen::MatrixXi jointE(numJoints, 2);
            for (unsigned int i = 0; i < numJoints; ++i)
            {
                const Eigen::Vector3f p0 = joints[i]->body0->q * joints[i]->r0 + joints[i]->body0->x;
                const Eigen::Vector3f p1 = joints[i]->body1->q * joints[i]->r1 + joints[i]->body1->x;

                jointP.row(2 * i) = p0;
                jointP.row(2 * i + 1) = p1;
                jointE.row(i) = Eigen::Vector2i(2 * i, 2 * i + 1);
            }

            auto pointCloud = polyscope::registerPointCloud("jointsPoint", jointP);
            pointCloud->setPointColor({ 0.0f, 0.0f, 1.0f });
            pointCloud->setPointRadius(0.005);
            auto curves = polyscope::registerCurveNetwork("jointsCurve", jointP, jointE);
            curves->setRadius(0.002f);
        }
    }

}

SimViewer::SimViewer() :
    m_dt(0.01f), m_subSteps(1), m_dynamicsTime(0.0f),
    m_paused(true), m_stepOnce(false),
    m_enableCollisions(true), m_enableScreenshots(false),
    m_resetState()
{
    m_resetState = std::make_unique<RigidBodySystemState>(*m_rigidBodySystem);
    reset();
}

SimViewer::~SimViewer()
{
}

void SimViewer::reset()
{
    std::cout << " ---- Reset ----- " << std::endl;
    m_resetState->restore(*m_rigidBodySystem);
    m_dynamicsTime = 0.0f;

    updateRigidBodyMeshes(*m_rigidBodySystem);
    polyscope::resetScreenshotIndex();
}

void SimViewer::save()
{
    std::cout << " ---- Saving current state ----- " << std::endl;
    m_resetState->save(*m_rigidBodySystem);
}

void SimViewer::start()
{
    // Setup Polyscope
    polyscope::options::programName = "MTI855 Devoir 02 Corps rigides";
    polyscope::options::verbosity = 0;
    polyscope::options::usePrefsFile = false;
    polyscope::options::alwaysRedraw = true;
    polyscope::options::ssaaFactor = 2;
    polyscope::options::openImGuiWindowForUserCallback = true;
    polyscope::options::groundPlaneHeightFactor = 0.0f; // adjust the plane height
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
    polyscope::options::buildGui = false;
    polyscope::options::maxFPS = -1;
    polyscope::options::groundPlaneEnabled = true;
    polyscope::options::screenshotExtension = ".png";

    // initialize
    polyscope::init();

    // Setup a viewing volume.
    polyscope::options::automaticallyComputeSceneExtents = false;
    polyscope::state::lengthScale = 10.0f;
    polyscope::state::boundingBox = std::tuple<glm::vec3, glm::vec3>{ {-5., 0, -5.}, {5., 5., 5.} };

    // Specify the update callback
    polyscope::state::userCallback = std::bind(&SimViewer::draw, this);

    // Add pre-step hook.
    m_rigidBodySystem->setPreStepFunc(std::bind(&SimViewer::preStep, this, std::placeholders::_1));

    // Show the window
    polyscope::show();

}

void SimViewer::drawGUI()
{
    ImGui::Text("Simulation:");
    ImGui::Checkbox("Pause", &m_paused);
    if (ImGui::Button("Step once"))
    {
        m_stepOnce = true;
    }
    if (ImGui::Button("Reset")) {
        reset();
    }
    if (ImGui::Button("Save")) {
        save();
    }

    ImGui::PushItemWidth(100);
    ImGui::SliderFloat("Time step", &m_dt, 0.0f, 0.1f, "%.3f");
    ImGui::SliderInt("Num. sub-steps", &m_subSteps, 1, 20, "%u");
    ImGui::SliderInt("Solver iters.", &(m_rigidBodySystem->solverIter), 1, 100, "%u");
    ImGui::SliderFloat("Friction coeff.", &(Contact::mu), 0.0f, 2.0f, "%.2f");
    ImGui::PopItemWidth();

    if (ImGui::Checkbox("Enable collision detecton", &m_enableCollisions)) {
        m_rigidBodySystem->setEnableCollisionDetection(m_enableCollisions);
    }

    ImGui::Checkbox("Enable screenshots", &m_enableScreenshots);

    if (ImGui::Button("Sphere-sphere distance")) {
        createSphereSphereDist();
    }
    if (ImGui::Button("Box on plane")) {
        createBoxOnPlane();
    }
    if (ImGui::Button("Sphere on box")) {
        createSphereOnBox();
    }
    if (ImGui::Button("Sphere inside box")) {
        createSphereInsideBox();
    }
    if (ImGui::Button("Sphere on cylinder 1")) {
        createCylinderOnSphere1();
    }
    if (ImGui::Button("Sphere on cylinder 2")) {
        createCylinderOnSphere2();
    }
    if (ImGui::Button("Rope ladder")) {
        createRopeLadder();
    }
    if (ImGui::Button("Marble box")) {
        createMarbleBox();
    }
    if (ImGui::Button("Custom scenario")) {
        createCustomScenario();
    }

    ImGui::Text("Step time: %3.3f ms", m_dynamicsTime);

}

void SimViewer::draw()
{
    drawGUI();

    updateJointViz(*m_rigidBodySystem);

    if( !m_paused || m_stepOnce )
    {
        auto start = std::chrono::high_resolution_clock::now();

        // Step the simulation.
        // The time step dt is divided by the number of sub-steps.
        //
        const float dt = m_dt / (float)m_subSteps;
        for(int i = 0; i < m_subSteps; ++i)
        {
            m_rigidBodySystem->step(dt);
        }
        auto stop = std::chrono::high_resolution_clock::now();

        updateRigidBodyMeshes(*m_rigidBodySystem);
        updateContactPoints(*m_rigidBodySystem);

        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        m_dynamicsTime = (float)duration.count() / 1000.0f;

        if (m_enableScreenshots)
        {
            polyscope::screenshot(false);
        }

        // Clear step-once flag.
        m_stepOnce = false;
    }
}

void SimViewer::createMarbleBox()
{
    Scenarios::createMarbleBox(*m_rigidBodySystem);
    m_resetState->save(*m_rigidBodySystem);
    updateRigidBodyMeshes(*m_rigidBodySystem);
    polyscope::resetScreenshotIndex();
}

void SimViewer::createSphereOnBox()
{
    Scenarios::createSphereOnBox(*m_rigidBodySystem);
    m_resetState->save(*m_rigidBodySystem);
    updateRigidBodyMeshes(*m_rigidBodySystem);
    polyscope::resetScreenshotIndex();
}

void SimViewer::createRopeLadder()
{
    Scenarios::createRopeLadder(*m_rigidBodySystem);
    m_resetState->save(*m_rigidBodySystem);
    updateRigidBodyMeshes(*m_rigidBodySystem);
    polyscope::resetScreenshotIndex();
}

void SimViewer::createCylinderOnSphere1()
{
    Scenarios::createCylinderSphereTest(*m_rigidBodySystem, Eigen::AngleAxisf::Identity());
    m_resetState->save(*m_rigidBodySystem);
    updateRigidBodyMeshes(*m_rigidBodySystem);
    polyscope::resetScreenshotIndex();
}

void SimViewer::createCylinderOnSphere2()
{
    Scenarios::createCylinderSphereTest(*m_rigidBodySystem, Eigen::AngleAxisf(-1.57f, Eigen::Vector3f(0.0f, 0.0f, 1.0f)));
    m_resetState->save(*m_rigidBodySystem);
    updateRigidBodyMeshes(*m_rigidBodySystem);
    polyscope::resetScreenshotIndex();
}

void SimViewer::createBoxOnPlane()
{
    Scenarios::createBoxOnPlane(*m_rigidBodySystem);
    m_resetState->save(*m_rigidBodySystem);
    updateRigidBodyMeshes(*m_rigidBodySystem);
    polyscope::resetScreenshotIndex();
}

void SimViewer::createSphereSphereDist()
{
    Scenarios::createSphereSphereDistance(*m_rigidBodySystem);
    m_resetState->save(*m_rigidBodySystem);
    updateRigidBodyMeshes(*m_rigidBodySystem);
    polyscope::resetScreenshotIndex();
}

void SimViewer::createSphereInsideBox()
{
    Scenarios::createSphereInsideBox(*m_rigidBodySystem);
    m_resetState->save(*m_rigidBodySystem);
    updateRigidBodyMeshes(*m_rigidBodySystem);
    polyscope::resetScreenshotIndex();
}

void SimViewer::createCustomScenario()
{
    Scenarios::createCustomScenario(*m_rigidBodySystem);
    m_resetState->save(*m_rigidBodySystem);
    updateRigidBodyMeshes(*m_rigidBodySystem);
    polyscope::resetScreenshotIndex();
}


void SimViewer::preStep(std::vector<RigidBody*>& _bodies)
{
    // do something useful here?
}
