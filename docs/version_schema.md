# Major.Schema.Results.Point format

## Outline
- Major: increment each time we seriously rework the pipeline
    - The v1 to v2 change qualified, though it’s not the minimum change that would qualify.
    - Expect to increment every 1-5 years.
- Schema: increment each time we restructure the outputs
    - The reworking with the addition of vertebrate infecting viruses would qualify
    - Consumers need to be updated because they will be looking for things in the wrong places
- Results: increment each time the results are no longer directly comparable with previous versions.
    - Historically a huge fraction of changes would require a bump here
    - So would lots of the changes we have planned
    - With the introduction of RBGD it starts to matter a lot that we have counts that can be compared across pipeline runs, and that we know when this is and isn’t the case.
- Point: increment for changes that don’t meet the criteria above
    - Changes that impact performance but not results, options that are off by default, and new outputs (it’s only changes to existing outputs that can’t go in Point releases)
    - Very small changes to outputs that don’t meaningfully impact RBGD or other pipeline consumers may also be considered for point releases, with team discussion.
    - I don’t think it is a problem if any individual segment (most likely Results) ends up very large.

## Notes
Only Major.Schema.Results goes into naming output directories, since consumers don’t need to consider Point (unless things have gone wrong).

We won’t want to reprocess all past data (which this time next year should be trillions of read pairs) each time we bump Results. So we introduce a concept of Stable, which is whatever version we currently are always running on new data and have run on all old data. We’ll need to be thoughtful about when to bump Stable, which will depend a lot on how expensive the pipeline ends up being to run. Sometimes it will make sense to backport bugfixes to Stable, and so Stable can have new Point releases. If a fix to Stable needed a new Results release, though, that would just be considering a new version to be Stable.

We would also introduce Unstable, which we run all new data through but don’t intend to (or haven’t decided if we will) run all past data through. This would be valuable for chimera detection and anything else that wants the best current analysis regardless of whether it’s consistent with history. I think this may just be Master, though, in which case we don’t need a new name.
