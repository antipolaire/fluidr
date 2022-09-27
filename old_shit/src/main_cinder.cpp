

#include "cinder/app/App.h"
#include "cinder/app/RendererGl.h"
#include "cinder/gl/gl.h"
#include "cinder/Rand.h"
#include "cinder/params/Params.h"

// #define SOLVER_JOS_STAM 1
#ifdef SOLVER_JOS_STAM
#include "solver_v2.hpp"
#else
#include "solver.hpp"
#endif

using namespace ci;
using namespace ci::app;
using namespace std;

static const int WINDOW_WIDTH = 600;
static const int WINDOW_HEIGHT = 600;

static const int FIELD_WIDTH = 200;
static const int FIELD_HEIGHT = 200;

static const float INITIAL_PARTICLE_SIZE = 2.0f;

struct Particle {
    vec3 home;
    vec3 pos; // position
    vec3 direction;
    float pressure;
};

class EvoCFP : public App {
public:
    // Allows us to override default window size, among other things.
    static void prepare(Settings *settings);

    void initInterface();

    void setup() override;

    void update() override;

    void draw() override;

    void cleanup() override;

    void keyDown(KeyEvent event) override;

    //! Called when the window is resized.
    void resize() override;

    void mouseDown(MouseEvent event) override;

    void mouseDrag(MouseEvent event) override;

private:

    vector<Particle> mParticles;

    gl::VboRef mParticleVbo;
    gl::BatchRef mParticleBatch;

#ifdef SOLVER_JOS_STAM
    FluidSolverJosStam *solver;
#else
    FluidSolver *solver;
#endif

    bool mPlay = false;

    float mFPS;

    vec2 mMouse = vec2(FIELD_WIDTH / 2., FIELD_WIDTH / 2.);

    params::InterfaceGlRef mParams;

    void initParticles();

    void updateParticles();

    gl::GlslProgRef createShaderProg(DataSourceRef, DataSourceRef, DataSourceRef);

};

void EvoCFP::keyDown(KeyEvent event) {
    switch (event.getCode()) {
        case KeyEvent::KEY_RIGHT:
            this->solver->fill(1);
            this->solver->step();
            break;
        case KeyEvent::KEY_SPACE:
            this->mPlay = !this->mPlay;
            break;
        case KeyEvent::KEY_F1:
            this->solver->fill(2);
            break;
        case KeyEvent::KEY_RETURN:
            break;
        case KeyEvent::KEY_F9:
            break;
        case KeyEvent::KEY_F10:
            break;
        case KeyEvent::KEY_F11:
            break;
        case KeyEvent::KEY_F12:
            break;
        default:
            break;
    }
}

void EvoCFP::prepare(Settings *settings) {
    settings->setWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
}

void EvoCFP::initInterface() {
    mParams = params::InterfaceGl::create("Parameters", ivec2(320, 100));
    mParams->setOptions("", "valueswidth=120");
    mParams->setOptions("", "refresh=0.05");
    mParams->addParam("FPS", &mFPS, false).step(0.1f);
    mParams->addParam("cos 0", &(this->solver->cos0), false).step(1.0f);
    mParams->addParam("cos 1", &(this->solver->cos1), false).step(0.1f);
    mParams->addParam("dt", &(this->solver->dt), false).step(0.001f);
#ifndef SOLVER_JOS_STAM
    mParams->addParam("dx", &(this->solver->dx), false).step(0.001f);
    mParams->addParam("rk_h", &(this->solver->rk_h), false).step(0.01f);
#endif
    mParams->addParam("viscosity", &(this->solver->viscosity), false).step(0.00001f);
    mParams->addParam("advect amount", &(this->solver->advectAmount), false).step(1.0f);
    mParams->maximize();
}

void EvoCFP::mouseDown(MouseEvent event) {
    mMouse.x = (float) event.getPos().x;
    mMouse.y = (float) event.getPos().y;
}

void EvoCFP::mouseDrag(MouseEvent event) {
    mMouse.x = (float) event.getPos().x;
    mMouse.y = (float) event.getPos().y;
}

void EvoCFP::cleanup() {
    delete this->solver;
}

void EvoCFP::setup() {

#ifdef SOLVER_JOS_STAM
    this->solver = new FluidSolverJosStam(FIELD_WIDTH, FIELD_HEIGHT);
#else
    this->solver = new FluidSolver(FIELD_WIDTH, FIELD_HEIGHT);
#endif

    initParticles();
    initInterface();

    // Allow the application to run at the same
    // frame rate as our monitor.
    gl::enableVerticalSync(true);
    disableFrameRate();
}

void EvoCFP::update() {
    mFPS = getAverageFps();

    // Use a fixed time step for a steady 60 updates per second,
    // regardless of our frame rate.
    static const double timestep = 1.0 / 60.0;

    // Keep track of time.
    static double time = getElapsedSeconds();
    static double accumulator = 0.0;

    // Calculate elapsed time since last frame.
    double elapsed = getElapsedSeconds() - time;
    time += elapsed;

    // Update stuff.
    accumulator += math<double>::min(elapsed, 0.1); // prevents 'spiral of death'
    while (accumulator >= timestep) {
        if(this->mPlay){
            this->solver->fill(1);
        }
        this->solver->step();
        updateParticles();
        accumulator -= timestep;
    }
}

void EvoCFP::draw() {
    gl::clear(ColorA::zero());
    gl::setMatricesWindowPersp(getWindowSize(), 60.0f, 1.0f, 10000.0f);
    gl::disableDepthRead();
    gl::disableDepthWrite();
    mParticleBatch->draw();
    mParams->draw();
}

void EvoCFP::resize() {
    initParticles();
}

void EvoCFP::initParticles() {
    mParticles.assign(FIELD_WIDTH * FIELD_HEIGHT, Particle());
    vec3 center = vec3(0.0);
    vec3 scale = vec3(getWindowWidth() / ((float) FIELD_WIDTH), getWindowHeight() / ((float) FIELD_HEIGHT), 1.0);
    for (int x = 0; x < FIELD_WIDTH; x++) {
        for (int y = 0; y < FIELD_HEIGHT; y++) {
            float p_x = x;
            float p_y = y;
            float p_z = 0.0f;
            int i = x + FIELD_WIDTH * y;
            auto &p = mParticles.at(i);
            p.home = vec3(p_x, p_y, p_z);
            p.pos = center + p.home * scale;
            p.direction = vec3(0.0, 0.0, 0.0);
            p.pressure = 0.0f;
        }
    }

    // Create particle buffer on GPU and copy over data.
    // Mark as streaming, since we will copy new data every frame.
    mParticleVbo = gl::Vbo::create(GL_ARRAY_BUFFER, mParticles, GL_STREAM_DRAW);

    // Describe particle semantics for GPU.
    geom::BufferLayout particleLayout;
    particleLayout.append(geom::Attrib::POSITION, 3, sizeof(Particle), offsetof(Particle, pos));
    particleLayout.append(geom::Attrib::CUSTOM_0, 3, sizeof(Particle), offsetof(Particle, direction));

    // Create mesh by pairing our particle layout with our particle Vbo.
    // A VboMesh is an array of layout + vbo pairs
    auto mesh = gl::VboMesh::create(mParticles.size(), GL_POINTS, {{particleLayout, mParticleVbo}});

    mParticleBatch = gl::Batch::create(mesh, createShaderProg(loadResource("passthr.vert"),
                                                              loadResource("passthr.frag"),
                                                              loadResource("dots.geom")),
                                       {{geom::CUSTOM_0, "direction"}});
    // gl::pointSize(INITIAL_PARTICLE_SIZE);

}

gl::GlslProgRef EvoCFP::createShaderProg(DataSourceRef vertexShaderFile, DataSourceRef fragmentShaderFile,
                                         DataSourceRef geometryShaderFile) {
    try {
        auto fmt = gl::GlslProg::Format()
                .vertex(vertexShaderFile)
                .fragment(fragmentShaderFile)
                .geometry(geometryShaderFile);

        return gl::GlslProg::create(fmt);
    }
    catch (const std::exception &exc) {
        console() << "Failed to load shader: " << exc.what() << std::endl;
    }
    return NULL;
}

void EvoCFP::updateParticles() {
    for (auto &p : mParticles) {
        p.direction.x = solver->getVelocityU(p.home.x, p.home.y);
        p.direction.y = solver->getVelocityV(p.home.x, p.home.y);
        p.direction.z = 0;
        p.pressure = solver->getPressure(p.home.x, p.home.y);
    }
    // Copy particle data onto the GPU.
    // Map the GPU memory and write over it.
    void *gpuMem = mParticleVbo->mapReplace();
    memcpy(gpuMem, mParticles.data(), mParticles.size() * sizeof(Particle));
    mParticleVbo->unmap();
}

CINDER_APP(EvoCFP, RendererGl(RendererGl::Options().msaa(8)), &EvoCFP::prepare)