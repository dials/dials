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

// The GitHub credential the branch source uses for scanning and status checks.
// It also posts the failure summary comment, so needs write access to pull
// requests (or issues) on dials/dials.
GITHUB_CREDENTIALS = 'cb3cd928-8897-4f7e-90f0-12d202e3c688'

// Posts a summary of the test failures to the pull request, as a single
// comment that later builds update. Nothing is posted for a passing build
// unless it needs to update an earlier failure comment. Jenkins is only
// reachable from the Diamond network, so this is how anyone else sees why a
// build failed. Written without backslashes, as Groovy would interpret them.
PR_COMMENT_SCRIPT = '''
import html
import json
import os
import sys
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET

MARKER = "<!-- jenkins-dials-regression -->"
NL = chr(10)
FENCE = "`" * 3
MAX_LISTED = 30
MAX_TRACEBACK = 3000
MAX_BODY = 60000

API = "https://api.github.com/repos/" + os.environ["GITHUB_REPO"]


def api(method, url, body=None):
    data = None if body is None else json.dumps(body).encode()
    request = urllib.request.Request(url, data=data, method=method)
    request.add_header("Authorization", "Bearer " + os.environ["GITHUB_TOKEN"])
    request.add_header("Accept", "application/vnd.github+json")
    with urllib.request.urlopen(request) as response:
        return json.load(response)


def find_failures():
    if not os.path.exists("output.xml"):
        return None
    failures = []
    for case in ET.parse("output.xml").iter("testcase"):
        for problem in case:
            if problem.tag in ("failure", "error"):
                name = case.get("classname", "") + "::" + case.get("name", "")
                failures.append((name, problem.text or problem.get("message") or ""))
    return failures


def make_body(status, failures):
    number = os.environ["BUILD_NUMBER"]
    lines = [MARKER]
    if status == "SUCCESS":
        lines.append("### :white_check_mark: Jenkins full regression passed")
    else:
        lines += ["### :x: Jenkins full regression: " + status.lower(), ""]
        if failures is None:
            lines.append("The build failed before any tests ran.")
        elif not failures:
            lines.append("No tests failed, so the build failed outside of the tests.")
        else:
            lines += ["{} test(s) failed:".format(len(failures)), ""]
            for i, (name, traceback) in enumerate(failures):
                if i == MAX_LISTED or sum(map(len, lines)) > MAX_BODY:
                    lines.append("...and {} more.".format(len(failures) - i))
                    break
                traceback = traceback[-MAX_TRACEBACK:].replace(FENCE, "'" * 3)
                lines += [
                    "<details><summary><code>{}</code></summary>".format(html.escape(name)),
                    "",
                    FENCE,
                    traceback.rstrip(),
                    FENCE,
                    "</details>",
                ]
    lines += [
        "",
        "Build [#{}]({}) (only reachable from the Diamond network)".format(
            number, os.environ["BUILD_URL"]
        ),
    ]
    return NL.join(lines)


def find_existing_comment(pr):
    existing = None
    page = 1
    while True:
        comments = api(
            "GET", API + "/issues/{}/comments?per_page=100&page={}".format(pr, page)
        )
        for comment in comments:
            if comment["body"].startswith(MARKER):
                existing = comment
        if len(comments) < 100:
            return existing
        page += 1


def main():
    pr = os.environ["CHANGE_ID"]
    status = os.environ["BUILD_STATUS"]
    body = make_body(status, find_failures())
    if os.environ.get("DRY_RUN"):
        print(body)
        return
    existing = find_existing_comment(pr)
    if existing:
        try:
            api("PATCH", existing["url"], {"body": body})
            return
        except urllib.error.HTTPError as e:
            # Someone else wrote a comment starting with the marker
            print("Could not update {}: {}".format(existing["html_url"], e))
    if status != "SUCCESS":
        api("POST", API + "/issues/{}/comments".format(pr), {"body": body})


main()
'''

def commentOnPullRequest() {
    writeFile file: 'pr_comment.py', text: PR_COMMENT_SCRIPT
    withCredentials([usernamePassword(credentialsId: GITHUB_CREDENTIALS, usernameVariable: 'GITHUB_USER', passwordVariable: 'GITHUB_TOKEN')]) {
        withEnv(["BUILD_STATUS=${currentBuild.currentResult}", 'GITHUB_REPO=dials/dials']) {
            def status = sh returnStatus: true, script: '''
                PYTHON=conda_base/bin/python
                [ -x "$PYTHON" ] || PYTHON=python3
                "$PYTHON" pr_comment.py
            '''
            if (status != 0) {
                echo "Could not comment on pull request #${env.CHANGE_ID}"
            }
        }
    }
}

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
                        conda_base/bin/pytest -v -rsxX -n auto --durations=0 --regression --junit-xml=output.xml ${WORKSPACE}/modules/dials
                        '''
                    }
                }
            }
            post {
                always {
                    junit   checksName: 'Tests',
                            testResults: 'output.xml',
                            allowEmptyResults: true
                    script {
                        if (env.CHANGE_ID) {
                            commentOnPullRequest()
                        }
                    }
                }
            }
        }
    }
}
