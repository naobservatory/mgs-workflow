#!/usr/bin/env python

#=============================================================================
# Imports
#=============================================================================

import os
import re
import sys
import argparse
from pathlib import Path
from typing import Dict, Set, List
from dataclasses import dataclass, field

#=============================================================================
# Dataclasses
#=============================================================================

@dataclass
class Process:
    name: str
    file_path: Path
    line_number: int
    module_name: str = None

@dataclass
class Module:
    name: str
    path: Path
    processes: Dict[str, Process] = field(default_factory = dict)

@dataclass
class Workflow:
    name: str
    file_path: Path
    is_anonymous: bool = False

#=============================================================================
# Main class
#=============================================================================

class NextflowAnalyzer:

    def __init__(self, pipeline_dir: str):
        # Directories
        self.pipeline_dir = Path(pipeline_dir)
        self.workflow_dir = self.pipeline_dir / "workflows"
        self.subworkflow_dir = self.pipeline_dir / "subworkflows" / "local"
        self.module_dir = self.pipeline_dir / "modules" / "local"
        # Initialize dictionaries
        self.workflows: Dict[str, Workflow] = {}
        self.modules: Dict[str, Module] = {}
        self.standalone_processes: Dict[str, Process] = {}
        self.workflow_dependencies: Dict[str, Dict[str, Set[str]]] = {}
        self.unused_components: Dict[str, Set[str]] = {}
        # Scan contents and construct lists
        self._scan_files()
        self._analyze_dependencies()
        self._find_unused_components()

    def _scan_files(self):
        """Scan directory for Nextflow files and identify workflows, modules and processes."""
        # First check for main.nf in root directory
        main_nf = self.pipeline_dir / "main.nf"
        if main_nf.exists():
            self.workflows["main"] = Workflow(name="main", file_path=main_nf, is_anonymous=True)
        # Scan workflows directory
        if self.workflow_dir.exists():
            for workflow_file in self.workflow_dir.rglob("*.nf"):
                with open(workflow_file) as f:
                    content = f.read()
                if re.search(r'workflow\s*{', content): # No anonymous workflows except main
                    raise ValueError(f"Unexpected anonymous workflow in file: {workflow_file}")
                workflow_matches = re.finditer(r'workflow\s+(\w+)\s*{', content)
                for match in workflow_matches:
                    workflow_name = match.group(1)
                    self.workflows[workflow_name] = Workflow(name=workflow_name, file_path = workflow_file)
        # Scan subworkflows directory
        if self.subworkflow_dir.exists():
            for subworkflow_path in self.subworkflow_dir.glob("**/main.nf"):
                subworkflow_name = subworkflow_path.parent.name
                self.workflows[subworkflow_name] = Workflow(name=subworkflow_name, file_path = subworkflow_path)
        # Scan modules directory
        if self.module_dir.exists():
            for module_path in self.module_dir.glob("**/main.nf"):
                module_name = module_path.parent.name
                self.modules[module_name] = Module(name=module_name, path=module_path)
                self._scan_file_for_processes(module_path, module_name)
        # Scan for standalone processes (excluding modules and workflows directories)
        for file_path in self.pipeline_dir.rglob("*.nf"):
            if ("modules/local" not in str(file_path) and
                "workflows" not in str(file_path) and
                file_path != main_nf):
                self._scan_file_for_processes(file_path)

    def _scan_file_for_processes(self, file_path: Path, module_name: str = None):
        """Scan a file for process definitions."""
        with open(file_path) as f:
            content = f.readlines()
        for line_num, line in enumerate(content, 1):
            process_match = re.search(r'process\s+(\w+)\s*{', line)
            if process_match:
                process_name = process_match.group(1)
                process = Process(name=process_name, file_path=file_path,
                                  line_number=line_num, module_name=module_name)
                if module_name:
                    self.modules[module_name].processes[process_name] = process
                else:
                    self.standalone_processes[process_name] = process

    def _analyze_dependencies(self):
        """Extract dependencies from workflows in pipeline."""
        for workflow_name in self.workflows.keys():
            workflow = self.workflows[workflow_name]
            self.workflow_dependencies[workflow.name] = {
                    "modules": set(),
                    "processes": set(),
                    "workflows": set()
                    }
            with open(workflow.file_path) as f:
                content = f.read()
                depend_matches = re.finditer(r'\s*include\s+{\s*(.*?)\s*}\s+from\s+[\'"](.*)/(.*)[\'"]\s*', content)
                for match in depend_matches:
                    item = match.group(1)
                    dirpath = match.group(2)
                    parent = match.group(3)
                    if re.search("/workflows", dirpath):
                        self.workflow_dependencies[workflow_name]["workflows"].add(item)
                    elif re.search("/subworkflows/local", dirpath):
                        self.workflow_dependencies[workflow_name]["workflows"].add(parent)
                    elif re.search("/modules/local", dirpath):
                        self.workflow_dependencies[workflow_name]["modules"].add(parent)
                        self.workflow_dependencies[workflow_name]["processes"].add(item)
                    else:
                        self.workflow_dependencies[workflow_name]["processes"].add(item)

    def _find_unused_components(self):
        """Identify modules, processes and workflows that aren't called within the main workflow."""
        used_modules, used_processes, used_workflows = self._extract_dependencies("main", set(), set(), set())
        unused_modules   = {name for name in self.modules.keys() if name not in used_modules}
        unused_workflows = {name for name in self.workflows.keys() if name not in used_workflows}
        unused_processes = {name for name in self.standalone_processes.keys() if name not in used_processes}
        for module_name in self.modules.keys():
            module = self.modules[module_name]
            unused_processes.update([name for name in module.processes.keys() if name not in used_processes])
        self.unused_components["modules"] = unused_modules
        self.unused_components["workflows"] = unused_workflows
        self.unused_components["processes"] = unused_processes

    def _extract_dependencies(self, workflow_name, modules, processes, workflows):
        workflows.add(workflow_name)
        processes.update(self.workflow_dependencies[workflow_name]["processes"])
        modules.update(self.workflow_dependencies[workflow_name]["modules"])
        dependent_workflows = self.workflow_dependencies[workflow_name]["workflows"]
        for dependent_workflow in dependent_workflows:
            modules, processes, workflows = self._extract_dependencies(dependent_workflow, modules, processes, workflows)
        return modules, processes, workflows

#=============================================================================
# Report generation function
#=============================================================================


def generate_nextflow_report(analyzer: NextflowAnalyzer, output_file: str):
    """Generate a detailed human-readable report from an Analyzer object."""
    with open(output_file, "w") as f:
        # Header
        f.write("=========================================\n")
        f.write("=== NEXTFLOW PIPELINE ANALYSIS REPORT ===\n")
        f.write("=========================================\n\n")
        f.write(f"Pipeline directory: {analyzer.pipeline_dir.absolute()}\n\n")
        # Summary statistics
        f.write("==========================\n")
        f.write("=== Summary Statistics ===\n")
        f.write("==========================\n\n")
        f.write(f"Total Workflows: {len(analyzer.workflows)}\n")
        f.write(f"Total Modules: {len(analyzer.modules)}\n")
        total_processes = len(analyzer.standalone_processes) + sum(
            len(m.processes) for m in analyzer.modules.values()
        )
        f.write(f"Total Processes: {total_processes}\n\n")
        f.write(f"Unused Workflows: {len(analyzer.unused_components['workflows'])}\n")
        f.write(f"Unused Modules: {len(analyzer.unused_components['modules'])}\n")
        f.write(f"Unused Processes: {len(analyzer.unused_components['processes'])}\n\n")
        # Workflow summary
        report_workflows(analyzer, f)

def report_workflows(analyzer: NextflowAnalyzer, output_stream):
    """Write workflow section of Nextflow report."""
    output_stream.write("=================\n")
    output_stream.write("=== Workflows ===\n")
    output_stream.write("=================\n\n")
    # Order workflows to start with main
    workflow_order = []
    if "main" in analyzer.workflows:
        workflow_order.append("main")
    workflow_order.extend(sorted(name for name in analyzer.workflows if name not in workflow_order))
    for name in workflow_order:
        # Print basic workflow information
        workflow = analyzer.workflows[name]
        output_stream.write(f"Workflow: {name}\n")
        output_stream.write(f"Location: {workflow.file_path.relative_to(analyzer.pipeline_dir)}\n")
        # Print information about dependencies
        if name in analyzer.workflow_dependencies:
            deps = analyzer.workflow_dependencies[name]
            dep_workflows = sorted(deps["workflows"])
            if dep_workflows:
                output_stream.write("Calls (sub)workflows:\n")
                for workflow in dep_workflows:
                    output_stream.write(f"\t- {workflow}\n")
            else:
                output_stream.write("Calls no (sub)workflows.\n")
            dep_modules = sorted(deps["modules"])
            if dep_modules:
                output_stream.write("Uses modules:\n")
                for module in dep_modules:
                    output_stream.write(f"\t- {module}\n")
            else:
                output_stream.write("Uses no modules.\n")
            dep_processes = sorted(deps["processes"])
            if dep_processes:
                output_stream.write("Uses processes:\n")
                for process in dep_processes:
                    output_stream.write(f"\t- {process}\n")
            else:
                output_stream.write("Uses no processes.\n")
        # Collate workflows that call this one
        called_by = [
                workflow for workflow in analyzer.workflow_dependencies.keys()
                if name in analyzer.workflow_dependencies[workflow]["workflows"]
                ]
        if called_by:
            output_stream.write("Called by workflows:\n")
            for workflow in sorted(called_by):
                output_stream.write(f"\t- {workflow}\n")
        elif name not in ["main"]:
            output_stream.write("Not called by any workflow.\n")
        output_stream.write("\n")

# TODO: Add missing sections to report: modules, processes, unused components


#class NextflowAnalyzer:
#
#    def generate_report(self, output_file: str):
#        """Generate a detailed human-readable report of the analysis."""
#
#        with open(output_file, 'w') as f:
#
#            for name in workflow_order:
#
#            # Modules section
#            f.write("\n=== Modules ===\n")
#            for name, module in sorted(self.modules.items()):
#                f.write(f"\nModule: {name}\n")
#                f.write(f"Location: {module.path.relative_to(self.pipeline_dir)}\n")
#
#                # List all processes in this module
#                if module.processes:
#                    f.write("Contains processes:\n")
#                    for proc_name, process in sorted(module.processes.items()):
#                        usage_status = "UNUSED" if proc_name in unused['module_processes'].get(name, []) else "used"
#                        f.write(f"  - {proc_name} ({usage_status})\n")
#                else:
#                    f.write("Contains no processes\n")
#
#                # List workflows that use this module
#                using_workflows = [
#                    workflow for workflow, deps in self.workflow_dependencies.items()
#                    if name in deps['modules']
#                ]
#                if using_workflows:
#                    f.write("Used by workflows:\n")
#                    for workflow in sorted(using_workflows):
#                        f.write(f"  - {workflow}\n")
#                else:
#                    f.write("Not used by any workflow\n")
#
#            # Standalone processes section
#            if self.standalone_processes:
#                f.write("\n=== Standalone Processes ===\n")
#                for name, process in sorted(self.standalone_processes.items()):
#                    f.write(f"\nProcess: {name}\n")
#                    f.write(f"Location: {process.file_path.relative_to(self.pipeline_dir)}:{process.line_number}\n")
#
#                    # List workflows that use this process
#                    using_workflows = [
#                        workflow for workflow, deps in self.workflow_dependencies.items()
#                        if name in deps['processes']
#                    ]
#                    if using_workflows:
#                        f.write("Used by workflows:\n")
#                        for workflow in sorted(using_workflows):
#                            f.write(f"  - {workflow}\n")
#                    else:
#                        f.write("Not used by any workflow\n")
#
#            # Unused components summary
#            f.write("\n=== Unused Components Summary ===\n")
#
#            if unused['modules']:
#                f.write("\nUnused modules:\n")
#                for module in sorted(unused['modules']):
#                    f.write(f"  - {module}\n")
#            else:
#                f.write("\nNo unused modules found\n")
#
#            if unused['processes']:
#                f.write("\nUnused standalone processes:\n")
#                for process in sorted(unused['processes']):
#                    f.write(f"  - {process}\n")
#            else:
#                f.write("\nNo unused standalone processes found\n")
#
#            unused_module_processes = {
#                module: procs for module, procs in unused['module_processes'].items()
#                if procs
#            }
#            if unused_module_processes:
#                f.write("\nUnused processes within modules:\n")
#                for module, processes in sorted(unused_module_processes.items()):
#                    f.write(f"  Module {module}:\n")
#                    for process in sorted(processes):
#                        f.write(f"    - {process}\n")
#            else:
#                f.write("\nNo unused module processes found\n")
#
#            if unused['subworkflows']:
#                f.write("\nUnused subworkflows:\n")
#                for workflow in sorted(unused['subworkflows']):
#                    f.write(f"  - {workflow}\n")
#            else:
#                f.write("\nNo unused subworkflows found\n")

#=============================================================================
# Run script
#=============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description='Analyze Nextflow pipeline to find unused workflows, modules, and processes.',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-d", "--pipeline_dir", default=os.getcwd(),
                        help="Path to the Nextflow pipeline directory (default: current working directory)")
    parser.add_argument("-o", "--output_path", default="workflow_report.txt",
                        help="Path to output file for detailed report (default: workflow_report.txt)")
    return parser.parse_args()

def main():
    # Parse arguments
    args = parse_args()
    pipeline_dir = args.pipeline_dir
    output_path = args.output_path
    # Check pipeline directory exists
    if not os.path.isdir(pipeline_dir):
        print(f"Error: Directory '{pipeline_dir}' does not exist", file=sys.stderr)
        sys.exit(1)
    # Initialize analyzer object
    analyzer = NextflowAnalyzer(pipeline_dir)
    # Create report
    generate_nextflow_report(analyzer, output_path)

if __name__ == "__main__":
    main()
