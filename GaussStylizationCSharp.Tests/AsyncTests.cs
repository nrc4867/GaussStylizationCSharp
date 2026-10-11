// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

using GaussStylization;

namespace GaussStylizationCSharp.Tests;

public class AsyncTests
{
    private static readonly double[] Cube = GaussStylizer.CubeDirections();

    [Fact]
    public async Task Async_MatchesSync_AndReportsProgress()
    {
        var (v, f) = TestMeshes.IcoSphere(2);
        var expected = GaussStylizer.Stylize(v, f, Cube, null, null, 3);
        var reports = new List<GaussStylizationProgress>();
        var progress = new SyncProgress(reports.Add);
        var actual = await GaussStylizer.StylizeAsync(v, f, Cube, null, null, 3, null, progress);
        Assert.Equal(expected, actual);
        Assert.Equal(new[] { 1, 2, 3 }, reports.Select(r => r.CompletedIterations));
        Assert.Equal(1.0, reports[^1].Fraction);
    }

    [Fact]
    public async Task PreCanceledToken_ReturnsCanceledTask()
    {
        var (v, f) = TestMeshes.IcoSphere(1);
        using var cts = new CancellationTokenSource();
        cts.Cancel();
        await Assert.ThrowsAnyAsync<OperationCanceledException>(
            () => GaussStylizer.StylizeAsync(v, f, Cube, null, null, 5, null, null, cts.Token));
    }

    [Fact]
    public async Task Cancel_DuringRun_StopsEarly()
    {
        var (v, f) = TestMeshes.IcoSphere(4);
        using var cts = new CancellationTokenSource();
        var progress = new SyncProgress(_ => cts.Cancel()); // cancel after the first iteration
        var task = GaussStylizer.StylizeAsync(v, f, Cube, null, null, 1000, null, progress, cts.Token);
        await Assert.ThrowsAnyAsync<OperationCanceledException>(() => task);
        Assert.True(task.IsCanceled);
    }

    [Fact]
    public void InvalidArguments_ThrowSynchronously()
    {
        var (v, f) = TestMeshes.IcoSphere(1);
        Assert.Throws<ArgumentOutOfRangeException>(() => { _ = GaussStylizer.StylizeAsync(v, f, Cube, null, null, -1); });
    }

    private sealed class SyncProgress : IProgress<GaussStylizationProgress>
    {
        private readonly Action<GaussStylizationProgress> _onReport;
        public SyncProgress(Action<GaussStylizationProgress> onReport) => _onReport = onReport;
        public void Report(GaussStylizationProgress value) => _onReport(value);
    }
}
