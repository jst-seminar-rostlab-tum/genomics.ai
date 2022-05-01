import React, { useState, useEffect } from "react";
import CircularProgress from "@mui/material/CircularProgress";
import styles from "./institutionTeamList.module.css";
import { getInstitutionTeams } from "shared/services/mock/teams";
import { getInstitution } from "shared/services/mock/institutions";
import InstitutionTeamCard from "components/institutions/InstitutionTeamCard";

function InstitutionTeamList({ onLeft, institution }) {
  const [teams, setTeams] = useState([]);
  const [teamsLoaded, setTeamsLoaded] = useState(false);

  useEffect(async () => {
    setTeams(await getInstitutionTeams(institution.id));
    setTeamsLoaded(true);
  }, []);

  if (!teamsLoaded) {
    return <CircularProgress />;
  }

  return (
    <div className={styles.content}>
      {teams.length === 0 ? "No teams." : ""}
      {teams.map((team) => (
        <div key={team.id}>
          <InstitutionTeamCard
            team={team}
            onLeft={onLeft}
            institution={institution}
          />
          <div className={styles.cardSpacing} />
        </div>
      ))}
    </div>
  );
}

export default InstitutionTeamList;
