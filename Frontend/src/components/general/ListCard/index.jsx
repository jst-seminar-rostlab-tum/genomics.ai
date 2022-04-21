import React from 'react';
import Card from '@mui/material/Card';
import Avatar from '@mui/material/Avatar';
import styles from './listCard.module.css';
import stringColor from 'shared/utils/stringColor';

function InstitutionCard({
  imageURL, title, description, nextToTitle, trailing,
}) {
  return (
    <Card sx={{ border: '1px solid #c8c8c8', boxShadow: '0 0 4px rgba(111, 111, 111, 0.25)' }} elevation={0}>
      <div className={styles.cardContent}>
        <div className={styles.start}>
          <div className={styles.cardImageWrapper}>
            <Avatar
              src={imageURL}
              alt={title}
              sx={{ backgroundColor: stringColor(title || ''), width: 42, height: 42 }}
            >
              {(title || '?')[0]}
            </Avatar>
          </div>
          <div className={styles.text}>
            <div className={styles.titleRow}>
              <h3 className={styles.title}>{title}</h3>
              {nextToTitle || (<></>)}
            </div>
            <div className={styles.description}>
              {description && (
                <p>{description}</p>
              )}
            </div>
          </div>
        </div>
        {trailing && (
          <div className={styles.trailing}>
            {trailing}
          </div>
        )}
      </div>
    </Card>
  );
}

export default InstitutionCard;
