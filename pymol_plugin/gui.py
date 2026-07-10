"""FlexAID∆S PyMOL GUI Panel.

Qt-based interface for docking result visualization and analysis.
"""

from pathlib import Path

try:
    from pymol.Qt import QtWidgets
    from pymol import cmd
except ImportError:
    raise ImportError("PyMOL Qt bindings not available")

from . import results_adapter as ro_adapter
from .session import SESSION


class FlexAIDSPanel(QtWidgets.QDialog):
    """Main GUI panel for FlexAID∆S plugin."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("FlexAID∆S: Entropy-Driven Docking")
        self.setMinimumWidth(420)
        self.setMinimumHeight(560)

        self._setup_ui()
        self._connect_signals()

        self._mode_ids: list[int] = []

    def _setup_ui(self):
        """Construct GUI layout."""
        layout = QtWidgets.QVBoxLayout(self)

        file_group = QtWidgets.QGroupBox("Load Docking Results")
        file_layout = QtWidgets.QHBoxLayout()

        self.file_path_edit = QtWidgets.QLineEdit()
        self.file_path_edit.setPlaceholderText("Select FlexAID output directory...")

        self.browse_btn = QtWidgets.QPushButton("Browse")
        self.load_btn = QtWidgets.QPushButton("Load")
        self.load_btn.setEnabled(False)

        file_layout.addWidget(self.file_path_edit)
        file_layout.addWidget(self.browse_btn)
        file_layout.addWidget(self.load_btn)
        file_group.setLayout(file_layout)
        layout.addWidget(file_group)

        adapter_info = QtWidgets.QLabel(
            "Loads result directories through the flexaidds Python API "
            "(ensemble thermodynamic ledger + CF scoring proxy)."
        )
        adapter_info.setWordWrap(True)
        layout.addWidget(adapter_info)

        self.status_label = QtWidgets.QLabel("No results loaded.")
        self.status_label.setWordWrap(True)
        layout.addWidget(self.status_label)

        mode_group = QtWidgets.QGroupBox("Binding Modes")
        mode_layout = QtWidgets.QVBoxLayout()

        self.mode_list = QtWidgets.QListWidget()
        self.mode_list.setSelectionMode(QtWidgets.QListWidget.SingleSelection)
        mode_layout.addWidget(self.mode_list)

        mode_group.setLayout(mode_layout)
        layout.addWidget(mode_group)

        thermo_group = QtWidgets.QGroupBox("Thermodynamics")
        thermo_layout = QtWidgets.QFormLayout()

        self.free_energy_label = QtWidgets.QLabel("-")
        self.enthalpy_label = QtWidgets.QLabel("-")
        self.entropy_label = QtWidgets.QLabel("-")
        self.entropy_term_label = QtWidgets.QLabel("-")
        self.n_poses_label = QtWidgets.QLabel("-")

        thermo_layout.addRow("ΔG (kcal/mol):", self.free_energy_label)
        thermo_layout.addRow("ΔH (kcal/mol):", self.enthalpy_label)
        thermo_layout.addRow("S (kcal/mol·K):", self.entropy_label)
        thermo_layout.addRow("TΔS (kcal/mol):", self.entropy_term_label)
        thermo_layout.addRow("# Poses:", self.n_poses_label)

        thermo_group.setLayout(thermo_layout)
        layout.addWidget(thermo_group)

        viz_group = QtWidgets.QGroupBox("Visualization")
        viz_layout = QtWidgets.QVBoxLayout()

        self.show_ensemble_btn = QtWidgets.QPushButton("Show Pose Ensemble")
        self.show_ensemble_btn.setEnabled(False)

        self.color_cf_btn = QtWidgets.QPushButton("Color by CF")
        self.color_cf_btn.setEnabled(False)

        self.color_free_energy_btn = QtWidgets.QPushButton("Color by Free Energy")
        self.color_free_energy_btn.setEnabled(False)

        self.show_representative_btn = QtWidgets.QPushButton("Show Representative Only")
        self.show_representative_btn.setEnabled(False)

        self.print_details_btn = QtWidgets.QPushButton("Print Mode Details")
        self.print_details_btn.setEnabled(False)

        self.entropy_heatmap_btn = QtWidgets.QPushButton("Show Entropy Heatmap")
        self.entropy_heatmap_btn.setEnabled(False)

        self.animate_btn = QtWidgets.QPushButton("Animate Mode Transition")
        self.animate_btn.setEnabled(False)

        self.itc_plot_btn = QtWidgets.QPushButton("ITC Compensation Plot")
        self.itc_plot_btn.setEnabled(False)

        self.unload_btn = QtWidgets.QPushButton("Unload Results")
        self.unload_btn.setEnabled(False)

        viz_layout.addWidget(self.show_ensemble_btn)
        viz_layout.addWidget(self.color_cf_btn)
        viz_layout.addWidget(self.color_free_energy_btn)
        viz_layout.addWidget(self.show_representative_btn)
        viz_layout.addWidget(self.print_details_btn)
        viz_layout.addWidget(self.entropy_heatmap_btn)
        viz_layout.addWidget(self.animate_btn)
        viz_layout.addWidget(self.itc_plot_btn)
        viz_layout.addWidget(self.unload_btn)

        viz_group.setLayout(viz_layout)
        layout.addWidget(viz_group)

        nrg_group = QtWidgets.QGroupBox("Export")
        nrg_layout = QtWidgets.QVBoxLayout()

        self.launch_nrgsuite_btn = QtWidgets.QPushButton("Launch NRGSuite")
        self.export_to_nrg_btn = QtWidgets.QPushButton("Export Mode Table")
        self.export_to_nrg_btn.setEnabled(False)

        nrg_layout.addWidget(self.launch_nrgsuite_btn)
        nrg_layout.addWidget(self.export_to_nrg_btn)

        nrg_group.setLayout(nrg_layout)
        layout.addWidget(nrg_group)

        close_btn = QtWidgets.QPushButton("Close")
        close_btn.clicked.connect(self.close)
        layout.addWidget(close_btn)

    def _connect_signals(self):
        """Wire up button click handlers."""
        self.browse_btn.clicked.connect(self._browse_directory)
        self.file_path_edit.textChanged.connect(
            lambda t: self.load_btn.setEnabled(bool(t.strip()))
        )
        self.load_btn.clicked.connect(self._load_results)
        self.mode_list.itemSelectionChanged.connect(self._on_mode_selected)
        self.show_ensemble_btn.clicked.connect(self._show_pose_ensemble)
        self.color_cf_btn.clicked.connect(self._color_by_cf)
        self.color_free_energy_btn.clicked.connect(self._color_by_free_energy)
        self.show_representative_btn.clicked.connect(self._show_representative)
        self.print_details_btn.clicked.connect(self._print_mode_details)
        self.entropy_heatmap_btn.clicked.connect(self._show_entropy_heatmap)
        self.animate_btn.clicked.connect(self._animate_modes)
        self.itc_plot_btn.clicked.connect(self._itc_compensation_plot)
        self.launch_nrgsuite_btn.clicked.connect(self._launch_nrgsuite)
        self.export_to_nrg_btn.clicked.connect(self._export_mode_table)
        self.unload_btn.clicked.connect(self._unload_results)

    def _browse_directory(self):
        """Open file dialog to select FlexAID output directory."""
        directory = QtWidgets.QFileDialog.getExistingDirectory(
            self, "Select FlexAID Output Directory"
        )
        if directory:
            self.file_path_edit.setText(directory)

    def _viz_buttons(self):
        return (
            self.show_ensemble_btn, self.color_cf_btn,
            self.color_free_energy_btn, self.show_representative_btn,
            self.print_details_btn, self.entropy_heatmap_btn,
            self.animate_btn, self.itc_plot_btn, self.export_to_nrg_btn,
            self.unload_btn,
        )

    def _disable_viz_buttons(self):
        """Disable all visualization buttons (called before load attempts)."""
        for btn in self._viz_buttons():
            btn.setEnabled(False)

    def _enable_viz_buttons(self, enabled: bool = True):
        for btn in self._viz_buttons():
            btn.setEnabled(enabled)

    def _load_results(self):
        """Load docking results from the selected output directory."""
        self._disable_viz_buttons()
        self.mode_list.clear()
        self._mode_ids.clear()
        output_dir = Path(self.file_path_edit.text().strip())
        if not output_dir.exists():
            self.status_label.setText("Directory not found.")
            QtWidgets.QMessageBox.warning(
                self, "Error", f"Directory not found: {output_dir}"
            )
            return

        try:
            ro_adapter.load_docking_results(str(output_dir))
        except Exception as exc:
            self.status_label.setText(f"Load failed: {exc}")
            QtWidgets.QMessageBox.warning(
                self,
                "Load Failed",
                f"Could not load results with the Python adapter:\n{exc}",
            )
            return

        result = SESSION.result if SESSION.result is not None else ro_adapter._loaded_result
        if result is None or not result.binding_modes:
            self.status_label.setText("No binding modes found.")
            QtWidgets.QMessageBox.warning(
                self,
                "No Results",
                f"No docking results could be parsed in:\n{output_dir}",
            )
            return

        sorted_modes = sorted(
            result.binding_modes,
            key=lambda mode: (
                mode.free_energy if mode.free_energy is not None else float("inf"),
                mode.rank,
                mode.mode_id,
            ),
        )

        for mode in sorted_modes:
            dg_str = (
                f"{mode.free_energy:.2f} kcal/mol"
                if mode.free_energy is not None
                else "N/A"
            )
            cf_str = f"{mode.best_cf:.2f}" if mode.best_cf is not None else "N/A"
            label = (
                f"mode{mode.mode_id}: F={dg_str}  CF={cf_str}  "
                f"({mode.n_poses} poses)"
            )
            self.mode_list.addItem(label)
            self._mode_ids.append(mode.mode_id)

        self._enable_viz_buttons(True)
        self.status_label.setText(
            f"Loaded {result.n_modes} mode(s) from {output_dir.name}."
        )

        if self.mode_list.count():
            self.mode_list.setCurrentRow(0)

        QtWidgets.QMessageBox.information(
            self,
            "Loaded",
            f"Loaded {result.n_modes} binding modes from {output_dir.name}.",
        )

    def _unload_results(self):
        """Unload results and clear the panel."""
        try:
            ro_adapter.unload_results(delete_objects=1)
        except Exception as exc:
            QtWidgets.QMessageBox.warning(self, "Unload Failed", str(exc))
            return
        self.mode_list.clear()
        self._mode_ids.clear()
        self._disable_viz_buttons()
        for label in (
            self.free_energy_label, self.enthalpy_label, self.entropy_label,
            self.entropy_term_label, self.n_poses_label,
        ):
            label.setText("-")
        self.status_label.setText("Results unloaded.")

    def _selected_mode_id(self, *, warn: bool = True) -> int | None:
        row = self.mode_list.currentRow()
        if row < 0 or row >= len(self._mode_ids):
            if warn:
                QtWidgets.QMessageBox.warning(
                    self, "No Mode Selected", "Please select a binding mode first."
                )
            return None
        return self._mode_ids[row]

    def _find_mode(self, mode_id: int):
        mode = SESSION.get_mode(mode_id)
        if mode is not None:
            return mode
        result = ro_adapter._loaded_result
        if result is None:
            return None
        for m in result.binding_modes:
            if m.mode_id == mode_id:
                return m
        return None

    def _on_mode_selected(self):
        """Update thermodynamics panel when a binding mode is selected."""
        mode_id = self._selected_mode_id(warn=False)
        if mode_id is None:
            return

        mode = self._find_mode(mode_id)
        if mode is None:
            return

        temperature = mode.temperature
        if temperature is None:
            if SESSION.result is not None and SESSION.result.temperature is not None:
                temperature = SESSION.result.temperature
            elif ro_adapter._loaded_result is not None:
                temperature = ro_adapter._loaded_result.temperature
        entropy_term = (
            (mode.entropy * temperature)
            if (mode.entropy is not None and temperature is not None)
            else None
        )

        def _fmt(val, fmt=".4f"):
            return f"{val:{fmt}}" if val is not None else "N/A"

        self.free_energy_label.setText(_fmt(mode.free_energy))
        self.enthalpy_label.setText(_fmt(mode.enthalpy))
        self.entropy_label.setText(_fmt(mode.entropy, ".6f"))
        self.entropy_term_label.setText(_fmt(entropy_term))
        self.n_poses_label.setText(str(mode.n_poses))

    def _show_pose_ensemble(self):
        """Render all poses in the selected binding mode."""
        mode_id = self._selected_mode_id()
        if mode_id is not None:
            ro_adapter.show_binding_mode(mode_id, show_all=1)

    def _color_by_cf(self):
        """Color poses by CF score."""
        mode_id = self._selected_mode_id()
        if mode_id is not None:
            ro_adapter.color_mode_by_score(mode_id, metric="cf")

    def _color_by_free_energy(self):
        """Color poses by free energy."""
        mode_id = self._selected_mode_id()
        if mode_id is not None:
            ro_adapter.color_mode_by_score(mode_id, metric="free_energy")

    def _show_representative(self):
        """Show only the representative pose."""
        mode_id = self._selected_mode_id()
        if mode_id is not None:
            ro_adapter.show_binding_mode(mode_id, show_all=0)

    def _print_mode_details(self):
        """Print mode thermodynamic details to the PyMOL console."""
        mode_id = self._selected_mode_id()
        if mode_id is not None:
            ro_adapter.show_mode_details(mode_id)

    def _show_entropy_heatmap(self):
        """Render entropy heatmap for the selected binding mode."""
        mode_id = self._selected_mode_id()
        if mode_id is not None:
            from .entropy_heatmap import render_entropy_heatmap
            render_entropy_heatmap(mode_id)

    def _animate_modes(self):
        """Animate transition between two binding modes."""
        result = SESSION.result if SESSION.result is not None else ro_adapter._loaded_result
        if result is None or len(result.binding_modes) < 2:
            QtWidgets.QMessageBox.warning(
                self, "Not Enough Modes",
                "At least two binding modes are needed for animation."
            )
            return

        mode_id = self._selected_mode_id()
        if mode_id is None:
            return

        # Find the next mode in the sorted list
        ids = self._mode_ids
        current_idx = ids.index(mode_id) if mode_id in ids else 0
        next_idx = (current_idx + 1) % len(ids)
        next_mode_id = ids[next_idx]

        from .mode_animation import animate_binding_modes
        animate_binding_modes(mode_id, next_mode_id)

    def _itc_compensation_plot(self):
        """Generate enthalpy-entropy compensation plot."""
        from .itc_comparison import plot_enthalpy_entropy_compensation
        plot_enthalpy_entropy_compensation()

    def _launch_nrgsuite(self):
        """Launch NRGSuite plugin (if installed)."""
        try:
            cmd.do("nrgsuite")
        except Exception as exc:
            QtWidgets.QMessageBox.warning(
                self,
                "NRGSuite Not Available",
                f"Could not launch NRGSuite: {exc}\n\n"
                "Install from: https://github.com/NRGlab/NRGsuite",
            )

    def _export_mode_table(self):
        """Export the current read-only mode table as TSV."""
        result = SESSION.result if SESSION.result is not None else ro_adapter._loaded_result
        if result is None:
            QtWidgets.QMessageBox.warning(
                self, "No Results", "Load a results directory first."
            )
            return

        output_dir = self.file_path_edit.text().strip()
        export_file, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Save Mode Table",
            str(Path(output_dir) / "flexaids_mode_table.tsv"),
            "Tab-Separated Files (*.tsv);;Text Files (*.txt);;All Files (*)",
        )
        if not export_file:
            return

        out_path = Path(export_file)
        with open(out_path, "w", encoding="utf-8") as fh:
            fh.write(
                "mode_id\trank\tn_poses\tbest_cf\tfree_energy\tenthalpy\tentropy\theat_capacity\ttemperature\n"
            )
            for mode in sorted(result.binding_modes, key=lambda m: (m.rank, m.mode_id)):
                fh.write(
                    f"{mode.mode_id}\t{mode.rank}\t{mode.n_poses}\t{mode.best_cf}\t"
                    f"{mode.free_energy}\t{mode.enthalpy}\t{mode.entropy}\t"
                    f"{mode.heat_capacity}\t{mode.temperature or result.temperature}\n"
                )

        QtWidgets.QMessageBox.information(
            self, "Exported", f"Mode table written to:\n{out_path}"
        )
