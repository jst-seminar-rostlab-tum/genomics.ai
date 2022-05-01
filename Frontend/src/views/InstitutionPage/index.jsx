import React, { useState, useEffect } from "react";
import { useParams } from "react-router-dom";
import TextField from "@mui/material/TextField";
import { getInstitution } from "shared/services/mock/institutions";
import InstitutionMemberList from "components/institutions/InstitutionMemberList";
import styles from "./institutionPage.module.css";
import Avatar from "@mui/material/Avatar";
import CircularProgress from "@mui/material/CircularProgress";
import { getInstitutionTeams } from "shared/services/mock/teams";
import InstitutionTeamList from "components/institutions/InstitutionTeamList";
import { useAuth } from "shared/context/authContext";

function InstitutionPage() {
  let { id } = useParams();
  id = parseInt(id, 10);
  const [institution, setInstitution] = useState({});
  const [user] = useAuth();
  const [institutionLoaded, setInstitutionLoaded] = useState(false);

  function isAdmin() {
    return (institution.adminIds || []).includes(user._id);
  }

  const handleDescriptionChange = (event) => {
    setInstitution({
      ...institution,
      description: event.target.value,
    });
  };

  useEffect(async () => {
    setInstitution(await getInstitution(id));
    setInstitutionLoaded(true);
  }, []);

  function onLeft(team) {
    setTeams(teams.filter((i) => i.id !== team.id));
  }

  if (!institutionLoaded) {
    return <CircularProgress />;
  }

  return (
    <>
      <div
        className={styles.background}
        style={{
          backgroundImage: `url(${institution.backgroundPictureURL})`,
          resizeMode: "stretch",
        }}
      >
        <div className={styles.institutionIcon}>
          <Avatar
            src={institution.avatarUrl}
            sx={{ width: 200, height: 200 }}
          />
        </div>
        <h1 className={styles.imageText}>
          <span>{institution.name}</span>
        </h1>
        <h3 className={styles.imageText}>
          <span>{institution.country}</span>
        </h3>
        <p className={styles.imageText}>
          <span>
            {institution.memberIds?.length + institution.adminIds?.length}
            {` Members`}
          </span>
        </p>
      </div>
      <div className={styles.test}>
        <section>
          <h2>Description</h2>
          <hr />
          <TextField
            id="description"
            multiline
            minRows={3}
            maxRows={5}
            value={institution.description}
            InputProps={{
              readOnly: !isAdmin(),
            }}
            style={{ width: "100%" }}
            onChange={handleDescriptionChange}
            variant="standard"
          />
        </section>
        <section>
          <h2>Teams</h2>
          <hr />
          <div className={styles.content}>
            <InstitutionTeamList
              onLeft={(t) => onLeft(t)}
              institution={institution}
            />
            <div className={styles.cardSpacing} />
          </div>
        </section>
        <section>
          <h2>Members</h2>
          <hr />
          <InstitutionMemberList
            institution={institution}
            // eslint-disable-next-line no-shadow
            onRemoved={(institution, removedMember) => {
              setInstitution({
                ...institution,
                adminIds: institution.adminIds.filter(
                  (mId) => mId !== removedMember.id
                ),
                memberIds: institution.memberIds.filter(
                  (mId) => mId !== removedMember.id
                ),
              });
            }}
          />
        </section>
      </div>
    </>
  );
}

export default InstitutionPage;
