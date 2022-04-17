import React, { useState, useEffect } from 'react';
import HeaderView from 'components/HeaderView';
import InstitutionCard from 'components/InstitutionCard';
import styles from './institutionOverview.module.css';
import queryMyInstitutions from 'shared/services/mock/institutions';

function InstitutionOverview({ sidebarShown }) {
  const [institutions, setInstitutions] = useState([]);
  useEffect(() => {
    queryMyInstitutions()
      .then((newInstitutions) => setInstitutions(newInstitutions))
      .catch((ignored) => { console.log(ignored); });
  }, [setInstitutions]);

  function onLeft(institution) {
    setInstitutions(institutions.filter((i) => i.id !== institution.id));
  }

  return (
    <HeaderView
      sidebarShown={sidebarShown}
      title="My Institutions"
    >
      <div className={styles.content}>
        {institutions.length === 0 ? 'No institutions.' : ''}
        {institutions.map((institution) => (
          <div key={institution.id}>
            <InstitutionCard
              institution={institution}
              onLeft={(inst) => onLeft(inst)}
            />
            <div className={styles.cardSpacing} />
          </div>
        ))}
      </div>
    </HeaderView>
  );
}

export default InstitutionOverview;
