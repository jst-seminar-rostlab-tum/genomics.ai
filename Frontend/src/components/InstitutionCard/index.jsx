import React, { useState, useEffect } from 'react';
import ListCard from 'components/ListCard';
import styles from './institutionCard.module.css';

import getUser from 'shared/services/mock/user';

function InstitutionCard({ institution }) {
  const [user, setUser] = useState({});
  useEffect(() => {
    getUser(institution.id)
      .then(setUser);
  }, [setUser]);

  const {
    name, profilePictureURL, adminIds,
  } = institution;
  return (
    <ListCard
      title={name}
      imageURL={profilePictureURL}
      nextToTitle={(
        <span className={styles.accessRightIndicator}>
          {adminIds.indexOf(user.id) !== -1 ? 'admin' : 'member'}
        </span>
      )}
    />
  );
}

export default InstitutionCard;
