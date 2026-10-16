def DATE = new Date().format("yyyy-MM-dd")

// Pull requests from forks only build for these GitHub users. This file is
// loaded by the Remote Jenkinsfile Provider from the dials/dials `jenkins`
// branch, so a pull request cannot change it (or supply its own Jenkinsfile).
ALLOWED_FORK_AUTHORS = [
    'aaronfinke',
    'dwmoreau',
    'jbeilstenedmands',
    'kmdalton',
    'ndevenish',
    'pblowey',
    'PrinceWalnut',
]

// CHANGE_FORK is only set for pull requests from forks. Pull requests from
// branches on dials/dials need write access to push, so are always allowed.
// CHANGE_AUTHOR comes from the GitHub API, so cannot be spoofed by the author.
boolean isAllowed() {
    if (!env.CHANGE_FORK) {
        return true
    }
    def author = env.CHANGE_AUTHOR?.toLowerCase()
    return author && ALLOWED_FORK_AUTHORS.any { it.toLowerCase() == author }
}

pipeline {
    // No agent until the author has been checked, so nothing from the pull
    // request is checked out or run beforehand.
    agent none

    options {
        quietPeriod 0
    }

    triggers {
        pollSCM ''
    }

    parameters {
        booleanParam(name: 'CLEAN_BUILD', defaultValue: false, description: 'Start from a fresh build. This will recompile DIALS from scratch.')
    }

    stages {
        stage('Check Author') {
            when {
                not { expression { isAllowed() } }
            }
            steps {
                echo "Not building: pull request from fork by ${env.CHANGE_AUTHOR}, who is not in the allowlist"
                script {
                    currentBuild.result = 'NOT_BUILT'
                }
            }
        }
        stage('Build and Test') {
            when {
                beforeAgent true
                expression { isAllowed() }
            }
            agent { label "dials && linux && rhel8" }
            stages {
                stage('Prepare') {
                    steps {
                        script {
                            if (!fileExists("${WORKSPACE}/${DATE}") || params.CLEAN_BUILD) {
                                deleteDir()
                                writeFile file: "${WORKSPACE}/${DATE}", text: ""
                            }
                        }
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
    }
}
