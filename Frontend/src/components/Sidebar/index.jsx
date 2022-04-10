import React, { useEffect, useState } from 'react';
import { Link } from 'react-router-dom';
import { CgFileDocument } from 'react-icons/all';
import dateFormat from 'dateformat';
import styles from './sidebar.module.css';
import arrowClose from '../../assets/arrow-close.png';
import geneIcon from '../../assets/gene.png';
import { filterInProgress, queryJobs } from '../../shared/services/StatusQueueLogic';
import { PROJECTS_UPDATE_INTERVAL } from '../../shared/utils/common/constants';

const Sidebar = ({ sidebarShown, toggleSidebar }) => {
  // COLLAPSED
  if (sidebarShown) {
    return (
      <div className={styles.sidebarNavCollapsed}>
        <div className={styles.sidebarWrap}>
          <div>
            <input
              className={styles.projectsBannerCollapsed}
              type="image"
              alt="toggle sidebar"
              src={geneIcon}
              onClick={toggleSidebar}
            />
          </div>
        </div>
      </div>
    );
  }

  const [jobs, setJobs] = useState([]);

  useEffect(() => {
    const updateJobs = () => queryJobs().then((newJobs) => setJobs(filterInProgress(newJobs)))
      .catch((ignored) => { console.log(ignored); });
    updateJobs();
    const intervalId = setInterval(updateJobs, PROJECTS_UPDATE_INTERVAL);
    return (() => clearInterval(intervalId));
  }, [setJobs]);

  // NOT COLLAPSED
  return (
    <div className={styles.sidebarNav}>
      <div className={styles.sidebarWrap}>
        <div style={{ flex: 1, flexDirection: 'row', alignItems: 'center' }}>
          <div style={{ paddingBlock: '35px', paddingLeft: '260px' }}>
            <input
              type="image"
              className={styles.toggleButton}
              alt="toggle-icon"
              src={arrowClose}
              onClick={toggleSidebar}
            />
          </div>
        </div>

        <div
          className={styles.projectsBanner}
        >
          <img
            alt="gene-icon"
            src={geneIcon}
            style={{ height: '35px', paddingLeft: '5px', paddingRight: '30px' }}
          />
          Projects
        </div>
        { jobs
          .filter((job) => job.status === 'DONE')
          .map((job) => (
            <Link
              className={styles.dropdownLink}
              to={{ pathname: '/result', search: `tsv=${encodeURIComponent(job.location)}` }}
              target="_blank"
              rel="noopener noreferrer"
              style={{ textDecoration: 'none' }}
              key={job.id}
            >
              <CgFileDocument style={{ color: 'white' }} />
              <span className={styles.sidebarLabel}>
                {dateFormat(new Date(job.uploadDate), 'dd/mm/yyyy hh:MM')}
              </span>
            </Link>
          ))}
      </div>
    </div>
  );
};

export default Sidebar;
