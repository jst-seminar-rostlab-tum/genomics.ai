import React, { useState, useEffect } from "react";
import CircularProgress from "@mui/material/CircularProgress";
import MemberList from "components/members/MemberList";
import InstitutionMemberRemoveButton from "components/institutions/InstitutionMemberRemoveButton";
import getUser from "shared/services/mock/user";
import styles from "./institutionMemberList.module.css";

function InstitutionTeamList({ team }) {
  const [user, setUser] = useState({});
  useEffect(() => {
    getUser().then(setUser);
  }, [setUser]);

  if (!institution.adminIds?.length || !institution.memberIds?.length) {
    return <CircularProgress />;
  }

  return (
    <MemberList
      memberIds={[...institution.adminIds, ...institution.memberIds]}
      nextToNameBuilder={(member) => (
        <span className={styles.accessRightIndicator}>
          {institution.adminIds.indexOf(member.id) !== -1 ? "admin" : "member"}
        </span>
      )}
      trailingBuilder={(member) =>
        institution.adminIds.includes(user.id) && user.id === member.id ? null : (
          <InstitutionMemberRemoveButton institution={institution} member={member} />
        )
      }
    />
  );
}

export default InstitutionTeamList;
