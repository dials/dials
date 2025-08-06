pipeline {
    agent { label "dials && linux && rhel8" }

    options {
        quietPeriod 0
    }

    triggers {
        pollSCM ''
    }

    stages {
        stage('Prepare') {
            steps {
                dir("modules") {
                    dir("dials") {
                        checkout scm
                    }
                }
                sh 'echo ccache >> modules/dials/.conda-envs/linux.txt'
                sh 'echo ninja >> modules/dials/.conda-envs/linux.txt'
                sh 'modules/dials/installer/bootstrap.py update base'
            }
        }
        stage('Build') {
            steps {
                // Always delete the build, and rely on ccache to make the build fast
                dir("build") {
                    deleteDir()
                }
                sh '''
                export CMAKE_GENERATOR=Ninja
                export CMAKE_CXX_COMPILER_LAUNCHER=${WORKSPACE}/conda_base/bin/ccache
                modules/dials/installer/bootstrap.py build
                '''
            }
        }
        stage('Test') {
            steps {
                sh '''
                export PATH="${WORKSPACE}/conda_base/bin:$PATH"
                conda_base/bin/pytest -v -rsxX -n auto --durations=0 --junit-xml=output.xml ${WORKSPACE}/modules/dials
                '''
            }
        }
    }

    post {
        always {
            junit   checksName: 'Tests',
                    testResults: 'output.xml'
        }
    }
}