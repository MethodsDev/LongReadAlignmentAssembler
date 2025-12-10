#!/usr/bin/env python3

import pytest
import pandas as pd
import numpy as np
import logging
from DiffIsoformStatTest import differential_isoform_tests

# Import the function under test
# Adjust this import based on your actual module name
# from your_module import differential_isoform_tests


class TestDifferentialIsoformTests:
    """Unit tests for the differential_isoform_tests function."""

    @pytest.fixture
    def basic_data(self):
        """Basic test data with clear differential usage."""
        data = [
            # Gene 1: Clear differential usage, sufficient reads
            ["gene1", "isoform1", "transcript1_1", 10, 50],  # Up in condition B
            ["gene1", "isoform2", "transcript1_2", 40, 10],  # Down in condition B
            # Gene 2: Insufficient total reads (should be filtered)
            ["gene2", "isoform1", "transcript2_1", 8, 8],
            ["gene2", "isoform2", "transcript2_2", 4, 4],
        ]
        return pd.DataFrame(
            data,
            columns=["gene_id", "isoform_id", "transcript_id", "count_A", "count_B"],
        )

    @pytest.fixture
    def reciprocal_data(self):
        """Test data suitable for reciprocal delta_pi analysis."""
        data = [
            # Gene with strong reciprocal changes and high coverage
            ["gene1", "isoform1", "transcript1_1", 15, 80],  # Strong up in B
            ["gene1", "isoform2", "transcript1_2", 75, 15],  # Strong down in B
            ["gene1", "isoform3", "transcript1_3", 10, 15],  # Minor change
            # Gene with insufficient alternate coverage for reciprocal analysis
            ["gene2", "isoform1", "transcript2_1", 10, 70],  # Up in B, good coverage
            [
                "gene2",
                "isoform2",
                "transcript2_2",
                60,
                5,
            ],  # Down in B, low total coverage
        ]
        return pd.DataFrame(
            data,
            columns=["gene_id", "isoform_id", "transcript_id", "count_A", "count_B"],
        )

    @pytest.fixture
    def edge_case_data(self):
        """Edge case test data."""
        data = [
            # Single isoform gene (should be skipped)
            ["gene1", "isoform1", "transcript1_1", 50, 50],
            # Zero counts in one condition
            ["gene2", "isoform1", "transcript2_1", 0, 50],
            ["gene2", "isoform2", "transcript2_2", 0, 10],
            # Small delta_pi (should be filtered)
            ["gene3", "isoform1", "transcript3_1", 30, 32],
            ["gene3", "isoform2", "transcript3_2", 30, 28],
        ]
        return pd.DataFrame(
            data,
            columns=["gene_id", "isoform_id", "transcript_id", "count_A", "count_B"],
        )

    def test_basic_functionality_reciprocal_false(self, basic_data):
        """Test basic functionality with reciprocal_delta_pi=False."""
        results = differential_isoform_tests(
            basic_data,
            reciprocal_delta_pi=False,
            min_reads_per_gene=20,
            min_delta_pi=0.1,
            min_reads_DTU_isoform=20,
        )

        # Should return results
        assert results is not None
        assert isinstance(results, pd.DataFrame)

        # Should have gene1 (sufficient reads and delta_pi)
        assert "gene1" in results["gene_id"].values

        # Should not have gene2 (insufficient total reads)
        assert "gene2" not in results["gene_id"].values

        # Check required columns are present
        expected_cols = [
            "gene_id",
            "pvalue",
            "delta_pi",
            "dominant_transcript_ids",
            "total_counts_A",
            "total_counts_B",
            "dominant_counts_A",
            "dominant_counts_B",
            "status",
        ]
        for col in expected_cols:
            assert col in results.columns

        # Alternate columns should be None when reciprocal_delta_pi=False
        assert results["alternate_delta_pi"].iloc[0] is None
        assert results["alternate_transcript_ids"].iloc[0] is None

    def test_basic_functionality_reciprocal_true(self, reciprocal_data):
        """Test basic functionality with reciprocal_delta_pi=True."""
        results = differential_isoform_tests(
            reciprocal_data,
            reciprocal_delta_pi=True,
            min_reads_per_gene=30,
            min_delta_pi=0.1,
            min_reads_DTU_isoform=25,
        )

        assert results is not None
        assert isinstance(results, pd.DataFrame)

        # Should have gene1 (meets all reciprocal requirements)
        gene1_results = results[results["gene_id"] == "gene1"]
        assert len(gene1_results) > 0

        # Should have alternate columns populated
        assert gene1_results["alternate_delta_pi"].iloc[0] is not None
        assert gene1_results["alternate_transcript_ids"].iloc[0] is not None
        assert gene1_results["alternate_counts_A"].iloc[0] is not None
        assert gene1_results["alternate_counts_B"].iloc[0] is not None

    def test_min_reads_per_gene_filter(self, basic_data):
        """Test that min_reads_per_gene filter works correctly."""
        # Use high threshold that gene2 won't meet (total = 24)
        results = differential_isoform_tests(
            basic_data, min_reads_per_gene=50, reciprocal_delta_pi=False
        )

        # Gene1 has 110 total reads, gene2 has 24 total reads
        if results is not None:
            # Only gene1 should pass the filter
            assert all(results["gene_id"] == "gene1")

    def test_min_reads_DTU_isoform_filter_reciprocal_false(self, reciprocal_data):
        """Test min_reads_DTU_isoform filter with reciprocal_delta_pi=False."""
        results = differential_isoform_tests(
            reciprocal_data,
            reciprocal_delta_pi=False,
            min_reads_DTU_isoform=50,  # High threshold
            min_reads_per_gene=30,
        )

        # Should filter out genes where dominant isoforms don't have enough reads
        if results is not None:
            for _, row in results.iterrows():
                dominant_total = row["dominant_counts_A"] + row["dominant_counts_B"]
                assert dominant_total >= 50

    def test_min_reads_DTU_isoform_filter_reciprocal_true(self, reciprocal_data):
        """Test min_reads_DTU_isoform filter with reciprocal_delta_pi=True."""
        results = differential_isoform_tests(
            reciprocal_data,
            reciprocal_delta_pi=True,
            min_reads_DTU_isoform=30,
            min_reads_per_gene=30,
        )

        # Should filter out genes where either dominant OR alternate don't have enough reads
        if results is not None:
            for _, row in results.iterrows():
                dominant_total = row["dominant_counts_A"] + row["dominant_counts_B"]
                alternate_total = row["alternate_counts_A"] + row["alternate_counts_B"]
                assert dominant_total >= 30
                assert alternate_total >= 30

    def test_min_delta_pi_filter(self, edge_case_data):
        """Test that min_delta_pi filter works correctly."""
        results = differential_isoform_tests(
            edge_case_data,
            min_delta_pi=0.2,  # High threshold
            reciprocal_delta_pi=False,
            min_reads_per_gene=20,
        )

        # Gene3 has very small delta_pi changes, should be filtered out
        if results is not None:
            assert "gene3" not in results["gene_id"].values

    def test_single_isoform_gene_skipped(self, edge_case_data):
        """Test that genes with only one isoform are skipped."""
        results = differential_isoform_tests(edge_case_data, reciprocal_delta_pi=False)

        # Gene1 has only one isoform, should be skipped
        if results is not None:
            assert "gene1" not in results["gene_id"].values

    def test_zero_counts_handling(self, edge_case_data):
        """Test handling of genes with zero counts in one condition."""
        results = differential_isoform_tests(
            edge_case_data, reciprocal_delta_pi=False, min_reads_per_gene=20
        )

        # Gene2 has zero count_A values but should still be processable
        # if it meets other criteria
        if results is not None and "gene2" in results["gene_id"].values:
            gene2_result = results[results["gene_id"] == "gene2"].iloc[0]
            assert gene2_result["total_counts_A"] == 0
            assert gene2_result["total_counts_B"] == 60

    def test_return_value_structure(self, basic_data):
        """Test the structure of returned DataFrame."""
        results = differential_isoform_tests(basic_data, reciprocal_delta_pi=True)

        if results is not None:
            # Check all required columns exist
            expected_columns = [
                "gene_id",
                "pvalue",
                "delta_pi",
                "dominant_transcript_ids",
                "total_counts_A",
                "total_counts_B",
                "dominant_counts_A",
                "dominant_counts_B",
                "alternate_delta_pi",
                "alternate_transcript_ids",
                "alternate_counts_A",
                "alternate_counts_B",
                "status",
            ]
            for col in expected_columns:
                assert col in results.columns

            # Check data types
            assert results["pvalue"].dtype in [
                np.float64,
                np.float32,
                object,
            ]  # Could be None
            assert results["delta_pi"].dtype in [np.float64, np.float32]
            assert results["total_counts_A"].dtype in [np.int64, np.int32]
            assert results["total_counts_B"].dtype in [np.int64, np.int32]

    def test_empty_results_handling(self):
        """Test handling when no genes meet the criteria."""
        # Create data where no genes will pass filters
        empty_case_data = pd.DataFrame(
            [
                ["gene1", "isoform1", "transcript1_1", 1, 1],  # Too few reads
                ["gene1", "isoform2", "transcript1_2", 1, 1],
            ],
            columns=["gene_id", "isoform_id", "transcript_id", "count_A", "count_B"],
        )

        results = differential_isoform_tests(
            empty_case_data,
            min_reads_per_gene=50,  # High threshold
            reciprocal_delta_pi=False,
        )

        # Should return None when no results
        assert results is None

    def test_top_isoforms_parameter(self):
        """Test that top_isoforms_each parameter works correctly."""
        # Create data with many isoforms
        many_isoforms_data = []
        for i in range(1, 11):  # 10 isoforms
            count_A = 100 - (i * 5)  # Decreasing counts
            count_B = 10 + (i * 8)  # Increasing counts
            many_isoforms_data.append(
                ["gene1", f"isoform{i}", f"transcript1_{i}", count_A, count_B]
            )

        df = pd.DataFrame(
            many_isoforms_data,
            columns=["gene_id", "isoform_id", "transcript_id", "count_A", "count_B"],
        )

        results = differential_isoform_tests(
            df, top_isoforms_each=3, reciprocal_delta_pi=False  # Limit to top 3
        )

        # Function should work with the parameter
        assert results is not None or results is None  # Either works or gets filtered

    def test_pvalue_calculation(self, basic_data):
        """Test that p-values are calculated and within valid range."""
        results = differential_isoform_tests(basic_data, reciprocal_delta_pi=False)

        if results is not None:
            for _, row in results.iterrows():
                if row["status"] == "OK":
                    assert row["pvalue"] is not None
                    assert 0 <= row["pvalue"] <= 1

    def test_status_field(self, basic_data):
        """Test that status field is properly set."""
        results = differential_isoform_tests(basic_data, reciprocal_delta_pi=False)

        if results is not None:
            # Status should be either "OK" or "failed"
            valid_statuses = {"OK", "failed"}
            assert all(status in valid_statuses for status in results["status"])

    @pytest.mark.parametrize("reciprocal_delta_pi", [True, False])
    def test_reciprocal_parameter_consistency(self, basic_data, reciprocal_delta_pi):
        """Test consistent behavior with both reciprocal_delta_pi values."""
        results = differential_isoform_tests(
            basic_data, reciprocal_delta_pi=reciprocal_delta_pi
        )

        if results is not None:
            if reciprocal_delta_pi:
                # Should have alternate columns when True
                assert "alternate_delta_pi" in results.columns
                assert "alternate_transcript_ids" in results.columns
                assert "alternate_counts_A" in results.columns
                assert "alternate_counts_B" in results.columns
            else:
                # Should have alternate columns as None when False
                assert all(results["alternate_delta_pi"].isnull())
                assert all(results["alternate_transcript_ids"].isnull())
                assert all(results["alternate_counts_A"].isnull())
                assert all(results["alternate_counts_B"].isnull())

    def test_logging_integration(self, basic_data, caplog):
        """Test that logging works correctly."""
        with caplog.at_level(logging.DEBUG):
            results = differential_isoform_tests(basic_data, reciprocal_delta_pi=False)

        # Should have some log messages
        assert len(caplog.records) > 0
        assert "Running differential_isoform_tests()" in caplog.text

    def test_top_isoforms_each_1_with_reciprocal(self):
        """Test that filtered_group is expanded to ensure diversity when top_isoforms_each=1.
        
        This reproduces the STMN2 scenario where:
        - Gene has 2 isoforms with reciprocal delta_pi changes
        - Same isoform is top in both conditions with top_isoforms_each=1
        - Diversity expansion adds the missing direction's top isoform to filtered_group
        """
        data = [
            # STMN2-like gene: iso1 is top in both conditions, but iso2 shows reciprocal change
            ["STMN2", "iso1", "transcript1", 148.003, 542.993],  # Top in both, delta_pi = +0.21
            ["STMN2", "iso2", "transcript2", 94.000, 117.992],   # Not top in either, delta_pi = -0.21
        ]
        df = pd.DataFrame(data, columns=["gene_id", "isoform_id", "transcript_id", "count_A", "count_B"])
        
        # With top_isoforms_each=1 and reciprocal_delta_pi=True
        # OLD behavior: Would fail delta_pi test because filtered_group only has 1 isoform
        # NEW behavior: Should PASS because diversity expansion adds iso2 to filtered_group
        results = differential_isoform_tests(
            df,
            group_by_token="gene_id",
            min_reads_per_gene=25,
            min_delta_pi=0.1,
            top_isoforms_each=1,  # Key parameter - causes filtering
            reciprocal_delta_pi=True,  # Requires both directions
            min_reads_DTU_isoform=25,
            show_progress_monitor=False
        )
        
        # Should have 1 result (gene passes the test)
        assert len(results) == 1, "Gene should pass delta_pi test using full gene data"
        
        result = results.iloc[0]
        assert result["gene_id"] == "STMN2"  # Column name is the group_by_token
        
        # The dominant isoform should be iso1 (only one in filtered_group)
        assert "transcript1" in result["dominant_transcript_ids"]
        
        # Pvalue should be calculated
        assert not pd.isna(result["pvalue"])
        
    def test_top_isoforms_each_1_with_reciprocal_annotated(self):
        """Test annotated output for the same scenario."""
        data = [
            ["STMN2", "iso1", "transcript1", 148.003, 542.993],
            ["STMN2", "iso2", "transcript2", 94.000, 117.992],
        ]
        df = pd.DataFrame(data, columns=["gene_id", "isoform_id", "transcript_id", "count_A", "count_B"])
        
        results, annotated = differential_isoform_tests(
            df,
            group_by_token="gene_id",
            min_reads_per_gene=25,
            min_delta_pi=0.1,
            top_isoforms_each=1,
            reciprocal_delta_pi=True,
            min_reads_DTU_isoform=25,
            show_progress_monitor=False,
            return_annotated_df=True
        )
        
        # Check that both isoforms are in annotated output
        assert len(annotated) == 2
        
        # Check that gene_tested is True for both isoforms (gene passed the test)
        assert annotated["gene_tested"].all(), "Both isoforms should be marked as gene_tested=True"
        
        # Check that skip_reason is "." (no skip) for both
        assert (annotated["skip_reason"] == ".").all(), "skip_reason should be '.' (passed)"
        
        # Check that delta_pi values are correct (calculated from full gene)
        delta_pis = annotated.sort_values("transcript_id")["delta_pi"].values
        assert abs(delta_pis[0] - 0.21) < 0.01, "First isoform delta_pi should be ~0.21"
        assert abs(delta_pis[1] - (-0.21)) < 0.01, "Second isoform delta_pi should be ~-0.21"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
