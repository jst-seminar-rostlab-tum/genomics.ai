import React, { useState, useEffect } from 'react';
import CircularProgress from '@mui/material/CircularProgress';
import MemberList from 'components/members/MemberList';
import InstitutionMemberRemoveButton from 'components/institutions/InstitutionMemberRemoveButton';
import InstitutionService from 'shared/services/Institution.service';
import styles from './institutionMemberList.module.css';
import { useAuth } from 'shared/context/authContext';

function InstitutionMemberList({ institution, onRemoved }) {
  const [user] = useAuth();
  const [members, setMembers] = useState([]);
  const [isLoading, setIsLoading] = useState(true);
  useEffect(async () => {
    setIsLoading(true);
    setMembers(await InstitutionService.getInstitutionMembers(institution.id));
    setIsLoading(false);
  }, []);

  if (isLoading) {
    return <CircularProgress />;
  }

  return (
    <MemberList
      members={members}
      nextToNameBuilder={(member) => (
        <span className={styles.accessRightIndicator}>
          {institution.adminIds.indexOf(member.id) !== -1 ? 'admin' : 'member'}
        </span>
      )}
      trailingBuilder={(member) => (
        user._id === member.id ? null : (
          <InstitutionMemberRemoveButton
            institution={institution}
            member={member}
            onRemoved={onRemoved}
          />
        )
      )}
    />
  );
}

export default InstitutionMemberList;
