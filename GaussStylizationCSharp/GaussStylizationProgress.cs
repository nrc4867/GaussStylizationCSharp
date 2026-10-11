// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
// Part of GaussStylization, a C# port of gaussStylization
// (https://github.com/ugogon/gaussStylization). See NOTICE.md.

namespace GaussStylization
{
    /// <summary>Progress of an asynchronous stylization, reported after each completed outer iteration.</summary>
    public readonly struct GaussStylizationProgress
    {
        /// <summary>Creates a progress report.</summary>
        public GaussStylizationProgress(int completedIterations, int totalIterations)
        {
            CompletedIterations = completedIterations;
            TotalIterations = totalIterations;
        }

        /// <summary>Outer iterations finished so far.</summary>
        public int CompletedIterations { get; }

        /// <summary>Outer iterations requested.</summary>
        public int TotalIterations { get; }

        /// <summary>Fraction done, in [0, 1].</summary>
        public double Fraction => TotalIterations <= 0 ? 1.0 : (double)CompletedIterations / TotalIterations;
    }
}
