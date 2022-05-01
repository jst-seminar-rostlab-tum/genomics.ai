import React, { useState, useEffect } from 'react';
import CircularProgress from '@mui/material/CircularProgress';
import MemberList from 'components/members/MemberList';
import TeamMemberRemoveButton from '../TeamMemberRemoveButton';
import getUser from 'shared/services/mock/user';
import styles from './teamMemberList.module.css';

function TeamMemberList({ team }) {
  const [user, setUser] = useState({});
  useEffect(() => {
    getUser().then(setUser);
  }, [setUser]);

  if (!team.adminIds?.length || !team.memberIds?.length) {
    return <CircularProgress />;
  }

  return (
    <MemberList
      memberIds={[...team.adminIds, ...team.memberIds]}
      nextToNameBuilder={(member) => (
        <span className={styles.accessRightIndicator}>
          {team.adminIds.indexOf(member.id) !== -1 ? 'admin' : 'member'}
        </span>
      )}
      trailingBuilder={(member) => (
        team.adminIds.includes(user.id) && user.id !== member.id ? (
          <TeamMemberRemoveButton team={team} member={member} />
        ) : null
      )}
    />
  );
}

export default TeamMemberList;
