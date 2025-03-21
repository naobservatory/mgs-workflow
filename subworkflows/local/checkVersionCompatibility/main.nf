/**********************************************
| SUBWORKFLOW: VERSION COMPATIBILITY CHECKING |
**********************************************/

/* Takes in pipeline and index versions plus respective compatibilities
and raises an error if they are incompatible. */

/**********************
| AUXILIARY FUNCTIONS |
**********************/

def isVersionLess(version1, version2) {
    /* Compare two semantic versions and return a boolean stating
    whether the first is less (i.e. older) than the second. */
    // Remove terminal tags (anything after a hyphen)
    def cleanVersion1 = version1.tokenize("-")[0]
    def cleanVersion2 = version2.tokenize("-")[0]
    // Split components by periods
    def v1Components = cleanVersion1.tokenize(".")
    def v2Components = cleanVersion2.tokenize(".")
    // Convert to integers (and raise error if unable)
    def v1IntComponents
    def v2IntComponents
    try {
        v1IntComponents = v1Components.collect{ it.toInteger() }
    } catch (NumberFormatException e) {
        def msg1 = "Invalid version format: version 1 (${version1}) contains non-integer components."
        throw new IllegalArgumentException(msg1)
    }
    try {
        v2IntComponents = v2Components.collect{ it.toInteger() }
    } catch (NumberFormatException e) {
        def msg2 = "Invalid version format: version 2 (${version2}) contains non-integer components."
        throw new IllegalArgumentException(msg2)
    }
    // Get the shortest version length
    def minLength = Math.min(v1IntComponents.size(), v2IntComponents.size())
    // Compare common components
    for (int i = 0; i < minLength; i++) {
        if (v1IntComponents[i] < v2IntComponents[i]) return true
        if (v1IntComponents[i] > v2IntComponents[i]) return false
    }
    // If components are equal, longer version is newer
    return v1IntComponents.size() < v2IntComponents.size()
}

/***********
| WORKFLOW |
***********/

workflow CHECK_VERSION_COMPATIBILITY {
    take:
        pipeline_version_path
        index_version_path
        pipeline_min_index_version_path
        index_min_pipeline_version_path
    main:
        // Read files and extract versions
        def pipeline_version = file(pipeline_version_path).readLines().first()
        def index_version = file(index_version_path).readLines().first()
        def pipeline_min_index_version = file(pipeline_min_index_version_path).readLines().first()
        def index_min_pipeline_version = file(index_min_pipeline_version_path).readLines().first()
        // Check version compatibilities
        if (isVersionLess(pipeline_version, index_min_pipeline_version)) {
            def msg_a = "Pipeline version is older than index minimum: ${pipeline_version} < ${index_min_pipeline_version}"
            throw new Exception(msg_a)
        }
        if (isVersionLess(index_version, pipeline_min_index_version)) {
            def msg_b = "Index version is older than pipeline minimum: ${index_version} < ${pipeline_min_index_version}"
            throw new Exception(msg_b)
        }
}
