import React, { useEffect, useState } from 'react';
import CircularProgress from '@mui/material/CircularProgress';
import { queryTeamJobs } from 'shared/services/mock/jobs';
import JobCard from '../JobCard';

function JobList({ teamId, forPart }) { // forPart can be "geneMapper" or "geneCruncher"
  const [jobs, setJobs] = useState([]);
  const [loading, setLoading] = useState(true);

  useEffect(async () => {
    setLoading(true);
    setJobs(await queryTeamJobs(teamId, forPart));
    setLoading(false);
  }, [teamId, forPart]);

  function loadingWrapper(orElse) {
    return loading ? (
      <CircularProgress />
    ) : orElse;
  }

  return (
    loadingWrapper(
      jobs.length ? (
        jobs.map((job) => (
          <JobCard key={job.id} job={job} />
        ))
      ) : (
        <span>No jobs.</span>
      ),
    )
  );
}

export default JobList;
