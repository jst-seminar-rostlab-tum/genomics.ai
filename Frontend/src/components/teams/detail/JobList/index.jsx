import React, { useEffect, useState } from 'react';
import CircularProgress from '@mui/material/CircularProgress';
import { queryTeamJobs } from 'shared/services/mock/jobs';
import styles from './jobList.module.css';
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
          <div key={job.id}>
            <JobCard job={job} />
            <div className={styles.spaceBetweenCards} />
          </div>
        ))
      ) : (
        <span>No jobs.</span>
      ),
    )
  );
}

export default JobList;
